//Joseph Zieg March 4 2019

#include "msa_utils.h"
#include "simfold.h"
#include <iostream>
#include <string.h>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>

bool generateConsensusSequence(std::vector<std::string> seqList, char *seq) {
    unsigned int i, n, s, nSeq;
    n = seqList[0].length();

    if (n > 0) {
        // check seqList for consistency
        for (s = 1; s < seqList.size(); s++) {
            if (seqList[s].length() != n) {
                fprintf(stderr, "Error generating consensus sequence: "
                                "Length of aligned sequence #%d does not match length of first sequence.\n\n", s + 1);
                return false;
            }
        }
        nSeq = s -
               1; // assumes a CLUSTAL type file parse, excludes the last sequence, being the conserved base sequence.

        for (i = 0; i < n; i++) {
            int fm = 0;
            char c;
            std::map<char, int> freq = {{'A', 0},
                                        {'C', 0},
                                        {'G', 0},
                                        {'U', 0}};

            for (s = 0; s < nSeq; s++) {
                freq[seqList[s][i]]++;
            }
            for (auto base : freq) {// find most frequent base
                if (base.second > fm && base.first != '-') {
                    c = base.first;
                    fm = base.second;
                }
            }
            seq[i] = c;
        }
    }
    return true;
}

bool generateConservationStructure(std::vector<std::string> seqList, char *struc) {
    //TODO: Take input conserved base sequence, convert to '_..._' type input for simfold,
    //      get best pairings back as struc.
    std::string conserved = seqList[seqList.size() - 1];
    int start, end;
    start = conserved.find('*', 0);
    end = conserved.rfind('*');
    char seq[end - start], sim_struc[end - start];
    int i;
    for (i = start; i < end; i++) {
        // Convert conserved cols to simfold input, truncate input sequence to increase fold likelihood.
        switch (conserved[i]) {
            case ' ':
                seq[i - start] = '.';
                break;
            case '*':
                seq[i - start] = '_';
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }
    simfold(seq, sim_struc);
    // Fill in full structure information for HFold.
    for (i = 0; i < conserved.length(); i++) {
        if (i < start) {
            struc[i] = '.';
        } else if (i >= start and i < end) {
            struc[i] = sim_struc[i - start];
        } else if (i >= end) {
            struc[i] = '.';
        }
    }
    return true;
}

format detectMSAValidFormat(char *path) {
    std::ifstream fp;
    std::string header; // char array for buffer dequeue
    size_t found; // flag for header verification.

    fp.open(path);
    if(!fp){
        return MSA_UNKNOWN;
    }
    std::getline(fp, header);
    found = header.find(">");
    if (found != std::string::npos) {
        return MSA_FASTA;
    }
    found = header.find("CLUSTAL");
    if (found != std::string::npos) {
        return MSA_CLUSTAL;
    }
    return MSA_UNKNOWN;
}

std::vector<std::string> readMSASequences(char *path, format type) {
    std::ifstream fp;
    std::string temp; // char array for buffer dequeue
    std::vector<std::string> seqList;
    // Focusing on CLUSTAL style input for the time being.
    switch (type) {
//        case MSA_FASTA:
//            seqList = readFASTASequences(path);
//            break;
        case MSA_CLUSTAL:
            seqList = readCLUSTALSequences(path);
        default:
            printf("Something went wrong reading the file.\n");
            break;
    }
    return seqList;
}

std::vector<std::string> readFASTASequences(char *path) {
    std::ifstream fp;
    std::string temp, seq; //char array for buffer dequeue
    std::vector<std::string> seqList;

    fp.open(path);
    if (!fp) {
        printf("File not found\n");
        exit(1);
    }
    while (getline(fp, temp).good()) {
        // iterate through FASTA file, copy individual seq lines into single string,
        // add to arr as new sequence comes up.
        if (temp.empty() || temp[0] == '>') {
            if (!seq.empty()) {
                seqList.emplace_back(seq);
                seq.clear();
            }
        } else if (!temp.empty()) {
            seq += temp;
        }
    }
    if (!seq.empty()) {
        seqList.emplace_back(seq);
        seq.clear();
    }
    fp.close();

    return seqList;
}

std::vector<std::string> readCLUSTALSequences(char *path) {
    // iterate through Clustal file, aggregate individual sequences into single string. Last sequence aggregate is the
    // conserved bases of the alignment.
    int seq_start;
    std::ifstream fp;
    std::string temp, id, seq; //char array for buffer dequeue
    std::vector<std::string> seqList;
    std::map<std::string, std::string> idMap;
    fp.open(path);
    if (!fp) {
        printf("File not found\n");
        exit(1);
    }
    while (getline(fp, temp).good()) {
        if (temp.find("CLUSTAL") == std::string::npos && !temp.empty()) {
            if (temp.at(0) != ' ') {
                seq_start = temp.rfind(' ') + 1;
                id = temp.substr(0, temp.find(' ', 0));
                seq = temp.substr(seq_start, std::string::npos);
                idMap[id] += seq;
            } else {
                idMap["STRUC"] += temp.substr(seq_start, std::string::npos);
            }
        }
    }
    for (auto pair : idMap) // find most frequent base
        seqList.emplace_back(pair.second);
    return seqList;
}
