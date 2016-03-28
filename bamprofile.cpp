#include <map>
#include <vector>
#include <string>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/SamReadGroup.h"
#include "fastahack/Fasta.h"
#include <cmath>

using namespace BamTools;
using namespace std;

#define PHRED_MAX 1000


short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double qualityChar2LongDouble(char c) {
    return static_cast<long double>(c) - 33;
}

long double lnqualityChar2ShortInt(char c) {
    return log(static_cast<short>(c) - 33);
}

char qualityInt2Char(short i) {
    return static_cast<char>(i + 33);
}

long double ln2log10(long double prob) {
    return M_LOG10E * prob;
}

long double log102ln(long double prob) {
    return M_LN10 * prob;
}

long double phred2ln(int qual) {
    return M_LN10 * qual * -.1;
}

long double ln2phred(long double prob) {
    return -10 * M_LOG10E * prob;
}

long double phred2float(int qual) {
    return pow(10, qual * -.1);
}

long double float2phred(long double prob) {
    if (prob == 1)
        return PHRED_MAX;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
        return PHRED_MAX;
    else
        return p;
}

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [-b FILE]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help            this dialog" << endl
         << "    -f, --fasta-reference FILE  the reference sequence" << endl
         << "    -b, --bam FILE        use this BAM as input (multiple allowed)" << endl
         << "    -r, --region REGION   limit alignments to those in this region (chr:start..end)" << endl
         << "    -e, --even-region N   output a samtools-style region delimiting every block of N alignments" << endl
         << endl
         << "Generates reports on the rate of putative mutations or errors in the input alignment data." << endl
         << "Alignments are read from the specified files, or stdin if none are specified" << endl
         << endl
         << "author: Erik Garrison <erik.garrison@gmail.com>" << endl;

}

void setRegion(BamMultiReader& reader, string& regionStr) {

    // parse the region string
    if (!regionStr.empty()) {

        map<string, int> refLength;
        map<string, int> refID;

        int id = 0;
        RefVector references = reader.GetReferenceData();
        for (RefVector::iterator r = references.begin(); r != references.end(); ++r) {
            refLength[r->RefName] = r->RefLength;
            refID[r->RefName] = id++;
        }

        // parse the region string
        string startSeq;
        int startPos;
        int stopPos;

        size_t foundFirstColon = regionStr.find(":");

        // we only have a single string, use the whole sequence as the target
        if (foundFirstColon == string::npos) {
            startSeq = regionStr;
            startPos = 0;
            stopPos = -1;
        } else {
            startSeq = regionStr.substr(0, foundFirstColon);
            size_t foundRangeDots = regionStr.find("..", foundFirstColon);
            if (foundRangeDots == string::npos) {
                startPos = atoi(regionStr.substr(foundFirstColon + 1).c_str());
                // differ from bamtools in this regard, in that we process only
                // the specified position if a range isn't given
                stopPos = startPos + 1;
            } else {
                startPos = atoi(regionStr.substr(foundFirstColon + 1, foundRangeDots - foundRangeDots - 1).c_str());
                // if we have range dots specified, but no second number, read to the end of sequence
                if (foundRangeDots + 2 != regionStr.size()) {
                    stopPos = atoi(regionStr.substr(foundRangeDots + 2).c_str()); // end-exclusive, bed-format
                } else {
                    stopPos = refLength[startSeq];
                }
            }
        }

        if (stopPos == -1) {
            stopPos = refLength[startSeq];
        }

        int startSeqRefID = refID[startSeq];

        if (!reader.LocateIndexes()) {
            cerr << "region specified, but could not open load BAM index" << endl;
            exit(1);
        } else {
            reader.SetRegion(startSeqRefID, startPos, startSeqRefID, stopPos);
        }

    }

}

int main(int argc, char** argv) {

    vector<string> inputFilenames;

    string regionStr;

    string fastaFile;

    int32_t even_region = 0;

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"region", required_argument, 0, 'r'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"even-region", required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:r:f:e:",
                         long_options, &option_index);

        if (c == -1)
            break;
 
        switch (c) {

        case '?':
            printUsage(argc, argv);
            return 0;
            break;

        case 'h':
            printUsage(argc, argv);
            return 0;
            break;

        case 'b':
            inputFilenames.push_back(optarg);
            break;

        case 'r':
            regionStr = optarg;
            break;

        case 'f':
            fastaFile = optarg;
            break;

        case 'e':
            even_region = atoi(optarg);
            break;

        default:
                return 1;
                break;
        }
    }

    if (fastaFile.empty()) {
        cerr << "no FASTA reference specified" << endl;
        return 1;
    }

    if (inputFilenames.empty()) {
        cerr << "no input files specified" << endl;
        return 1;
    }

    BamMultiReader reader;
    if (!reader.Open(inputFilenames)) {
        cerr << "could not open input BAM files" << endl;
        return 1;
    }

    setRegion(reader, regionStr);

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    FastaReference fr;
    fr.open(fastaFile);

    long unsigned int alignedBases = 0;
    long unsigned int mismatchCount = 0;
    long unsigned int gapCount = 0;
    map<int, long unsigned int> mismatches;
    map<int, long unsigned int> gaps;

    long int lowestReferenceBase = 0;
    long unsigned int referenceBases = 0;
    int currentRefSeqID = -1;

    map<short, uint64_t> qual_hist;

    string region_start_seq;
    int32_t region_start_pos = 0;
    int32_t region_end_pos = 0;
    int64_t seen = 0;

    BamAlignment al;
    while (reader.GetNextAlignment(al)) {
        if (al.IsMapped()) {

            long unsigned int endpos = al.GetEndPosition();
            // this happens when we switch reference sequences
            if (currentRefSeqID != al.RefID) {
                //cout << "al.GetEndPosition() = " << endpos << "  lowestReferenceBase = " << lowestReferenceBase << "  reset" << endl;
                currentRefSeqID = al.RefID;
                referenceBases += lowestReferenceBase;
                lowestReferenceBase = 0;
                if (!region_start_seq.empty()) {
                    cout << region_start_seq << ":" << region_start_pos << "-" << region_end_pos << endl;
                }
                region_start_seq = referenceIDToName[al.RefID];
                region_start_pos = 0;
                seen = 0;
            } else if (endpos > lowestReferenceBase) {
                //cout << "al.GetEndPosition() = " << endpos << "  lowestReferenceBase = " << lowestReferenceBase << "  adding " << endpos - lowestReferenceBase << endl;
                referenceBases += endpos - lowestReferenceBase;
                lowestReferenceBase = endpos;
            }

            region_end_pos = al.GetEndPosition();

            if (even_region) {
                if (++seen > even_region) {
                    // write the region
                    cout << region_start_seq << ":" << region_start_pos << "-" << region_end_pos << endl;
                    region_start_pos = al.Position;
                    region_end_pos = al.GetEndPosition();
                    seen = 0;
                }
                continue;
            }

            // record the qualities
            for (string::iterator c = al.Qualities.begin(); c != al.Qualities.end(); ++c) {
                ++qual_hist[qualityChar2ShortInt(*c)];
            }
            
            //cout << al.Position << endl;
            //cout << al.AlignedBases << endl;
            string refsequence = fr.getSubSequence(referenceIDToName[al.RefID], al.Position, al.GetEndPosition() - (al.Position - 1));
            //cout << refsequence << endl;

            alignedBases += al.QueryBases.size();

            int rp = 0; int sp = 0;
            vector<CigarOp>::const_iterator cigarIter = al.CigarData.begin();
            vector<CigarOp>::const_iterator cigarEnd  = al.CigarData.end();
            for ( ; cigarIter != cigarEnd; ++cigarIter ) {
                unsigned int l = cigarIter->Length;
                char t = cigarIter->Type;

                if (t == 'M') { // match or mismatch

                    int firstMismatch = -1;

                    for (int i=0; i<l; i++) {

                        // extract aligned base
                        char b = al.QueryBases.at(rp);

                        // get reference allele
                        char sb = refsequence.at(sp);

                        // record mismatch if we have a mismatch here
                        if (firstMismatch >= 0) {
                            if (b == sb) {
                                // mismatch termination
                                // register multi-base mismatch
                                int length = rp - firstMismatch;
                                //string qualstr = rQual.substr(rp - length, length);
                                ++mismatches[length];
                                mismatchCount += length;
                                firstMismatch = -1;
                            } else {
                                // mismatch extension
                            }
                        } else {
                            if (b != sb) {
                                // mismatch state
                                firstMismatch = rp;
                            } else {
                                // match state
                            }
                        }

                        // update positions
                        ++sp;
                        ++rp;
                    }

                } else if (t == 'D') {
                    ++gaps[-l];
                    ++gapCount;
                    sp += l;
                } else if (t == 'I') {
                    ++gaps[l];
                    ++gapCount;
                    rp += l;
                } else if (t == 'S') {
                    rp += l;
                } else if (t == 'H') {
                } else if (t == 'N') {
                    sp += l;
                    rp += l;
                }
            }
        }
    }

    reader.Close();

    if (even_region) {
        cout << region_start_seq << ":" << region_start_pos << "-" << region_end_pos << endl;        
        return 0;
    }
    
    cout << "reference bases:\t" << referenceBases << endl;
    cout << "total aligned bases:\t" << alignedBases << endl;
    cout << "mean alignment depth:\t" << (long double) alignedBases / (long double) referenceBases << endl;
    cout << "total mismatched bases:\t" << mismatchCount << endl;
    cout << "total gap bases:\t" << gapCount << endl;
    cout << "mismatch rate per aligned bp:\t" << (long double) mismatchCount / (long double) alignedBases << endl;
    cout << "gap rate per aligned bp:\t" << (long double) gapCount / (long double) alignedBases << endl;
    cout << "mismatch + gap per aligned bp:\t" << (long double) ( gapCount + mismatchCount ) / (long double) alignedBases << endl;
    cout << endl;

    cout << "mismatch length distribution" << endl;
    cout << "length\tcount\trate (per aligned base)" << endl;
    for (map<int, long unsigned int>::iterator p = mismatches.begin(); p != mismatches.end(); ++p) {
        cout << p->first << "\t" << p->second << "\t" << (long double) p->second / (long double) alignedBases << endl;
    }
    cout << endl;

    cout << "gap length distribution:" << endl;
    cout << "length\tcount\trate (per aligned base)" << endl;
    for (map<int, long unsigned int>::iterator p = gaps.begin(); p != gaps.end(); ++p) {
        cout << p->first << "\t" << p->second << "\t" << (long double) p->second / (long double) alignedBases << endl;
    }

    cout << endl;
    cout << "quality histogram of alignmed reads" << endl;
    cout << "#qual\tcount" << endl;
    for (map<short, uint64_t>::const_iterator q = qual_hist.begin(); q != qual_hist.end(); ++q) {
        cout << q->first << "\t" << q->second << endl;
    }
    //cout << endl;

    return 0;

}
