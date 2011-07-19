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

using namespace BamTools;
using namespace std;

void printUsage(int argc, char** argv) {

    cerr << "usage: " << argv[0] << " [-b FILE]" << endl
         << endl
         << "options:" << endl
         << "    -h, --help         this dialog" << endl
         << "    -b, --bam FILE     use this BAM as input (multiple allowed)" << endl
         << "    -r, --region REGION  limit alignments to those in this region (chr:start..end)" << endl
         << endl
         << "Generates reports on the rate of putative mutations or errors in the input alignment data." << endl
         << "Alignments are read from the specified files, or stdin if none are specified" << endl
         << endl
         << "author: Erik Garrison <erik.garrison@bc.edu>" << endl;

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

    // parse command-line options
    int c;

    while (true) {

        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bam",  required_argument, 0, 'b'},
            {"region", required_argument, 0, 'r'},
            {"fasta-reference", required_argument, 0, 'f'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hb:r:f:",
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
    unsigned int currentRefSeqID = 0;

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
            } else if (endpos > lowestReferenceBase) {
                //cout << "al.GetEndPosition() = " << endpos << "  lowestReferenceBase = " << lowestReferenceBase << "  adding " << endpos - lowestReferenceBase << endl;
                referenceBases += endpos - lowestReferenceBase;
                lowestReferenceBase = endpos;
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
    //cout << endl;

    return 0;

}
