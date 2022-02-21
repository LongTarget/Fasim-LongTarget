#include <iostream>
#include <string.h>
#include <sstream>
#include <fstream>
#include "ssw_cpp.h"
#include "ssw.h"
#include "sim.h"
#define N 50
using std::string;
using std::cout;
using std::endl;
using std::ifstream;

struct mycutregion
{
	mycutregion() {};
	mycutregion(int ei, int ej) :start(ei), end(ej) {};
	int start;
	int end;
};

void getAlignment(StripedSmithWaterman::Alignment alignment,
	string ref_seq,
	string read_seq,
	string ref_seq_src,
	const int8_t* table,
	string &ref_align,
	string &read_align,
	string &ref_align_src);

void convertMyTriplex(StripedSmithWaterman::Alignment alignment,
	std::vector<struct triplex> &triplex_list,
	string ref_seq,
	string read_seq,
	string read_seq_src,
	const int8_t* table,
	long dnaStartPos,
	long rule,
	long strand,
	long Para,
	int penaltyT,
	int penaltyC,
	int ntMin,
	int ntMax);

void cutSequence(string& seq, vector<string>& seqsVec, vector<int>& seqsStartPos,
	int cutLength, int overlapLength, int &cut_num)
{
	unsigned int pos = 0;
	int tmpa = 0;
	int tmpb = 0;
	seqsVec.clear();
	seqsStartPos.clear();
	string cutSeq;
	while (pos < seq.size())
	{
		cutSeq = seq.substr(pos, cutLength);
		seqsVec.push_back(cutSeq);
		seqsStartPos.push_back(pos);
		pos += cutLength;
		pos -= overlapLength;
		tmpa++;
	}
	cut_num = tmpa;
}

bool compMyTriplexSingle(const triplex &a, const triplex &b)
{
	return a.score > b.score;
}

bool compMyTriplexMultiple(const triplex &a, const triplex &b)
{
	// need to sort myTriplexList based on score.
	if (a.stari == b.stari)
	{
		if (a.starj == b.starj)
		{
			return a.score > b.score;
		}
		else
		{
			return a.starj > b.starj;
		}
	}
	else
	{
		return a.starj > b.starj;
	}

}

bool compMyTriplexMultiple2(const triplex &a, const triplex &b)
{
	// need to sort myTriplexList based on score.
	if (a.endi == b.endi)
	{
		if (a.starj == b.starj)
		{
			return a.score > b.score;
		}
		else
		{
			return a.starj < b.starj;
		}
	}
	else
	{
		return a.starj < b.starj;
	}

}

bool sameMyTriplex(const triplex &a, const triplex &b)
{
	if (a.stari == b.stari && a.starj == b.starj && a.endi == b.endi &&
		a.endj == b.endj && a.score == b.score)
	{
		return true;
	}
	else if (b.stari >= a.stari && b.starj >= a.starj && b.endi <= a.endi &&
		b.endj <= a.endj && b.score < a.score)
	{
		return true;
	}
	else
	{
		return false;
	}

}

void fastSIM(string& strA, string& strB, string& strSrc,
	long dnaStartPos, long min_score, float parm_M,
	float parm_I, float parm_O, float parm_E,
	vector<struct triplex>& triplex_list,
	long strand, long Para, long rule,
	int ntMin, int ntMax, int penaltyT,
	int penaltyC)
{
	int32_t maskLen = 15;
	StripedSmithWaterman::Aligner aligner;
	// Declares a default filter
	StripedSmithWaterman::Filter filter;
	// Declares an alignment that stores the result
	StripedSmithWaterman::Alignment alignment;
	struct mycutregion tmpregion;
	vector<struct mycutregion> mycut_list;
	vector<struct triplex> myTriplexList;
	// Aligns the query to the ref
	//aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
	//std::vector<struct scoreInfo> finalScoreInfo;
	std::vector<struct StripedSmithWaterman::scoreInfo> finalScoreInfo;
	std::vector<struct StripedSmithWaterman::scoreInfo> smallRevScoreInfo;
	std::vector<struct StripedSmithWaterman::scoreInfo> finalRevScoreInfo;
	std::vector<struct StripedSmithWaterman::scoreInfo> testFinalScore;
	int i = 0, j = 0, seqBegin = 0, seqEnd = 0;
	bool beginFound = false;
	int cutLength = 300;
	int midPoint = 0;
	int tmpLoop = 0;
	const int8_t nt_table[128] = {
	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 0, 4, 1,	4, 4, 4, 2,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 4, 4, 4,	3, 0, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 0, 4, 1,	4, 4, 4, 2,	4, 4, 4, 4,	4, 4, 4, 4,
	4, 4, 4, 4,	3, 0, 4, 4,	4, 4, 4, 4,	4, 4, 4, 4
	};
	aligner.preAlign(strA.c_str(), strB.c_str(), strB.size(), filter, &alignment,
		maskLen, min_score, finalScoreInfo, parm_M,parm_I);
//	for (i = 0; i < finalScoreInfo.size(); i++)
//	{
//		fprintf(stderr, "%d is score %d is position (final)\n",
//			finalScoreInfo[i].score, finalScoreInfo[i].position);
//	}
//	string queryRev = strA;
	string smallSeq;
//	string minorSeq;
//	string smallSeqRev;
//	reverseSeq(queryRev);
//	string refRev = strB;
//	reverseSeq(refRev);
//	aligner.preAlign(queryRev.c_str(), refRev.c_str(), refRev.size(),
//		filter, &alignment, maskLen, min_score, testFinalScore, parm_M,parm_I);

//	for (i = testFinalScore.size() - 1; i > -1; i--)
//	{
//		fprintf(stderr, "%d is score %d is position (test)\n",
//			testFinalScore[i].score, strB.size() - testFinalScore[i].position - 1);
//	}
	for (i = 0; i < finalScoreInfo.size(); i++)
		//if(1) // DEBUG.
	{
		//i = 4; // DEBUG.
		// Here we can use testFinalScore to help locate the begin position.
		//beginFound = false;
		//for(j = 0; j < testFinalScore.size(); j++)
		float Iden = 0.6;
		int cutlength;
		while(Iden<=1){
            cutlength = (int)(finalScoreInfo[i].score+24)/(9*Iden-4)+1;
            smallSeq = strB.substr((finalScoreInfo[i].position - cutlength + 1 > 0 ? finalScoreInfo[i].position - cutlength + 1 : 0), cutlength);
            aligner.Align(strA.c_str(), smallSeq.c_str(), smallSeq.size(), filter,&alignment, maskLen);
            if(alignment.sw_score==finalScoreInfo[i].score) break;
            Iden += 0.1;
		}
		alignment.ref_begin = alignment.ref_begin + finalScoreInfo[i].position - cutlength + 1;
        alignment.ref_end = alignment.ref_end + finalScoreInfo[i].position - cutlength + 1;
        convertMyTriplex(alignment,
                    myTriplexList,
                    strA,
                    strB,
                    strSrc,
                    nt_table,
                    dnaStartPos,
                    rule,
                    strand,
                    Para,
                    penaltyT,
                    penaltyC,
                    ntMin,
                    ntMax);

//		for (j = 0; j < testFinalScore.size(); j++)
//		{
//			if (finalScoreInfo[i].score == testFinalScore[j].score)
//			{
//				// first check score, then check position, position
//				// must be within cutLength.
//				if (finalScoreInfo[i].position -
//					(strB.size() - testFinalScore[j].position-1) > 0
//					&& finalScoreInfo[i].position -
//					(strB.size() - testFinalScore[j].position-1) < cutLength)
//				{
//					// found.
//					seqBegin = strB.size() - testFinalScore[j].position - 1;
//					seqEnd = finalScoreInfo[i].position + 1;
//					// THis should be carefully considered.
//					tmpregion = mycutregion(seqBegin, seqEnd);
//					mycut_list.push_back(tmpregion);
//					//beginFound = true;
//					//break;
//				}
//			}
//		}
//
//		if (mycut_list.size()>0)
//		{
//			// we can precisely get the alignment.
//			for (j = 0; j < mycut_list.size(); j++) {
//				tmpregion = mycut_list[j];
//				seqBegin = tmpregion.start;
//				smallSeq = strB.substr(tmpregion.start, tmpregion.end - seqBegin + 1);
//				aligner.Align(strA.c_str(), smallSeq.c_str(), smallSeq.size(), filter,&alignment, maskLen);
//				if (finalScoreInfo[i].score != alignment.sw_score)
//				{
////					 a better alignment is reported, try to find sub-optimal
////					 alignment whose score is identical with finalScoreInfo[i].score.
//					midPoint = 0;
//					while (alignment.sw_score != finalScoreInfo[i].score) // is this OK?
//					{
//						//tmpLoop++;
//						//if(tmpLoop > 5)
//						//{
//
//						//	break; // DEBUG.
//						//}
//
//						//cout << midPoint << " is midPoint " << alignment.ref_begin <<
//						//	" is ref_begin " << alignment.ref_end << " is ref_end" <<
//						//	endl;
//						midPoint = (alignment.ref_end + alignment.ref_begin) / 2
//							+ midPoint;
//						// get smallSeq again, which should be minorSeq.
//						minorSeq = smallSeq.substr(midPoint, smallSeq.size());
//						// will this be OK???
//						aligner.Align(strA.c_str(), minorSeq.c_str(),
//							minorSeq.size(), filter, &alignment, maskLen);
//						if (alignment.sw_score < min_score)
//						{
//
//							// if this happens, it means that this alignment
//							// contains a better alignment, which can't be
//							// recognized by current method.
//							// 3585-3713 score: 122 ==> 3632-3696 score: 127.
//							// then 3585-3713 can't be recognized.
//							break;
//						}
//						//smallSeq = minorSeq;
//					}
//
//					if (alignment.sw_score > min_score)
//					{
////						cout << alignment.sw_score << " is sw_score " <<
////							finalScoreInfo[i].position << " is position " <<
////							alignment.sw_score_next_best << " is sw_score_next_best " <<
////							alignment.ref_begin + seqBegin + 1 <<
////							" is ref_begin " << alignment.ref_end + seqBegin + 1
////							<< " is ref_end (beginFoundCalibrate)"
////							<< endl;
//						alignment.ref_begin = alignment.ref_begin + seqBegin + midPoint;
//						alignment.ref_end = alignment.ref_end + seqBegin + midPoint;
//							convertMyTriplex(alignment,
//								myTriplexList,
//								strA,
//								strB,
//								strSrc,
//								nt_table,
//								dnaStartPos,
//								rule,
//								strand,
//								Para,
//								penaltyT,
//								penaltyC,
//								ntMin,
//								ntMax);
//
//					}
//                }
//
//				else
//				{
//					if (alignment.sw_score > min_score)
//					{
////						cout << alignment.sw_score << " is sw_score " <<
////							finalScoreInfo[i].position << " is position " <<
////							alignment.sw_score_next_best << " is sw_score_next_best " <<
////							alignment.ref_begin + seqBegin + 1 <<
////							" is ref_begin " << alignment.ref_end +
////							seqBegin + 1 << " is ref_end (beginFound)" << endl;
//						alignment.ref_begin = alignment.ref_begin + seqBegin;
//						alignment.ref_end = alignment.ref_end + seqBegin;
//							convertMyTriplex(alignment,
//								myTriplexList,
//								strA,
//								strB,
//								strSrc,
//								nt_table,
//								dnaStartPos,
//								rule,
//								strand,
//								Para,
//								penaltyT,
//								penaltyC,
//								ntMin,
//								ntMax);
//					}
//
//				}
//			}
//			//smallSeq = strB.substr(seqBegin, seqEnd - seqBegin + 1);
//			//cout << "beginFound: " << seqBegin << " is seqBegin " <<
//			//seqEnd << " is seqEnd" << endl;
//		}
//		else
//		{
////		    int cutlength = (int)(finalScoreInfo[i].score+80)/1.4+1;
//			smallSeq = strB.substr((finalScoreInfo[i].position - cutLength + 1 > 0 ? finalScoreInfo[i].position -  + 1 : 0), cutLength);
//			aligner.Align(strA.c_str(), smallSeq.c_str(), smallSeq.size(), filter,&alignment, maskLen);
//			if (finalScoreInfo[i].position != alignment.ref_end +
//				(finalScoreInfo[i].position - cutLength + 1 > 0 ? finalScoreInfo[i].position - cutLength + 1 : 0))
//			{
//				// means that better alignment is reported, which may be reported
//				// in other segments, need to get midPoint and get smallSeq again.
//
//				// TOTHINK: will this happen again?
//				// Actually this will happen again...
//				//tmpLoop = 0;
//				//while(finalScoreInfo[i].position != alignment.ref_end
//				//	+ finalScoreInfo[i].position - minorSeq.size() + 1)
//				midPoint = 0;
//				while (alignment.sw_score != finalScoreInfo[i].score) // is this OK?
//				{
//					//tmpLoop++;
//					//if(tmpLoop > 5)
//					//{
//
//					//	break; // DEBUG.
//					//}
//
////					cout << midPoint << " is midPoint " << alignment.ref_begin <<
////						" is ref_begin " << alignment.ref_end << " is ref_end" <<
////						endl;
//					midPoint = (alignment.ref_end + alignment.ref_begin) / 2
//						+ midPoint;
//					// get smallSeq again, which should be minorSeq.
//					minorSeq = smallSeq.substr(midPoint, smallSeq.size());
//					// will this be OK???
//
//					aligner.Align(strA.c_str(), minorSeq.c_str(),
//						minorSeq.size(), filter, &alignment, maskLen);
//					if (alignment.sw_score < min_score)
//					{
//
//						// if this happens, it means that this alignment
//						// contains a better alignment, which can't be
//						// recognized by current method.
//						// 3585-3713 score: 122 ==> 3632-3696 score: 127.
//						// then 3585-3713 can't be recognized.
//						break;
//					}
//					//smallSeq = minorSeq;
//				}
//
//
//
////				cout << alignment.sw_score << " is sw_score " <<
////					finalScoreInfo[i].position << " is position " <<
////					alignment.sw_score_next_best << " is sw_score_next_best " <<
////					alignment.ref_begin + finalScoreInfo[i].position -
////					cutlength  + 1 + 1 << " is ref_begin " <<
////					alignment.ref_end + finalScoreInfo[i].position -
////					cutlength + 1 + 1 << " is ref_end (minorSeq)" << endl;
//				alignment.ref_begin = alignment.ref_begin + finalScoreInfo[i].position -
//					cutLength  + 1;
//				alignment.ref_end = alignment.ref_end + finalScoreInfo[i].position -
//					cutLength  + 1;
//
//					convertMyTriplex(alignment,
//								myTriplexList,
//								strA,
//								strB,
//								strSrc,
//								nt_table,
//								dnaStartPos,
//								rule,
//								strand,
//								Para,
//								penaltyT,
//								penaltyC,
//								ntMin,
//								ntMax);
//
//			}
//			else
//			{
////				cout << alignment.sw_score << " is sw_score " <<
////					finalScoreInfo[i].position << " is position " <<
////					alignment.sw_score_next_best << " is sw_score_next_best " <<
////					alignment.ref_begin + finalScoreInfo[i].position -
////					cutlength + 1 + 1 << " is ref_begin " <<
////					alignment.ref_end + finalScoreInfo[i].position -
////					cutlength + 1 + 1 << " is ref_end (smallSeq)" << endl;
//				alignment.ref_begin = alignment.ref_begin + finalScoreInfo[i].position -
//					cutLength + 1;
//				alignment.ref_end = alignment.ref_end + finalScoreInfo[i].position -
//					cutLength + 1;
//                    convertMyTriplex(alignment,
//								myTriplexList,
//								strA,
//								strB,
//								strSrc,
//								nt_table,
//								dnaStartPos,
//								rule,
//								strand,
//								Para,
//								penaltyT,
//								penaltyC,
//								ntMin,
//								ntMax);
//
//
//			}
//			// to check smallRevScoreInfo.
//			//smallSeqRev = smallSeq;
//			//reverseSeq(smallSeqRev);
//			//aligner.preAlign(queryRev.c_str(), smallSeqRev.c_str(),
//			//	smallSeqRev.size(), filter, &alignment, maskLen, threshold,
//			//	smallRevScoreInfo);
//			//for(j = 0; j < smallRevScoreInfo.size(); j++)
//			//{
//			//	fprintf(stderr, "%d is score %d is position %d is finalScore %d is finalPosition (smallSeqRev)\n",
//			//		smallRevScoreInfo[j].score,
//			//		smallRevScoreInfo[j].position,
//			//		finalScoreInfo[i].score,
//			//		finalScoreInfo[i].position);
//			//}
//			// clear smallRevScoreInfo.
//			//std::vector<struct StripedSmithWaterman::scoreInfo>().swap(smallRevScoreInfo);
//		}
		mycut_list.clear();
	}
	std::sort(myTriplexList.begin(), myTriplexList.end(), compMyTriplexMultiple);
	// then unique it.
	myTriplexList.erase(std::unique(myTriplexList.begin(), myTriplexList.end(),
		sameMyTriplex), myTriplexList.end());

	// need to sort again based on endi.
	std::sort(myTriplexList.begin(), myTriplexList.end(), compMyTriplexMultiple2);
	// and erase sub-alignment based on endi.
	myTriplexList.erase(std::unique(myTriplexList.begin(), myTriplexList.end(),
		sameMyTriplex), myTriplexList.end());

	// then sort again, based on score.
	std::sort(myTriplexList.begin(), myTriplexList.end(), compMyTriplexSingle);
	for (i = 0; i < (myTriplexList.size() > N ? N : myTriplexList.size()); i++) {
		triplex_list.push_back(myTriplexList[i]);
	}
	//triplex_list.assign(myTriplexList.begin(), myTriplexList.begin() + (myTriplexList.size() > N ? N : myTriplexList.size()));
}


void convertMyTriplex(StripedSmithWaterman::Alignment alignment,
	std::vector<struct triplex> &triplex_list,
	string read_seq,
	string ref_seq,
	string ref_seq_src,
	const int8_t* table,
	long dnaStartPos,
	long rule,
	long strand,
	long Para,
	int penaltyT,
	int penaltyC,
	int ntMin,
	int ntMax)
{
	int i = 0, j = 0, nLoop = 0;
	int match = 0, mis_match = 0;
	float identity = 0.0;
	int nt = 0;
	float tri_score = 0.0;
	float hashvalue = 0.0, prescore = 0.0;
	string ref_align;
	string read_align;
	string ref_align_src;
	int refStart, refEnd;
	//int strand = 0, Para = 1; //DEBUG.
	//int penaltyT = -1000, penaltyC = 0; // DEBUG.
	//int ntMin = 20, ntMax = 100000; // DEBUG.
	float score = 0.0;
	// main for loop to convert myTriplexList to SIM's triplex_list.
	// firstly need to get read_align and ref_align.
	getAlignment(alignment, ref_seq,  read_seq,  ref_seq_src, table, ref_align, read_align,ref_align_src);

//	std::cout << ref_align << " is ref_align" << std::endl;
//	std::cout << read_align << " is read_align" << std::endl;
	// then get identity.

	nt = ref_align.length();
	for (i = 0; i < nt; i++)
	{
		if (ref_align[i] == read_align[i])
		{
			match++;
		}
		else
		{
			mis_match++;
		}
	}
	identity = (float)(100 * match) / (float)(match + mis_match);
//	cout << identity << " is identity" << endl;
//	cout << ref_align << "is ref_align" << endl;
//	cout << read_align << "is read_align" << endl;
//	cout << ref_align_src << "is ref_align_src" << endl;
	// then calculate tri_score, will this be different with triplex_list
	// due to the different alignment?
	tri_score = 0.0;
	char prechar = 0, curchar = 0;
	//DEBUG, manually set ref_align and read_align.
	//ref_align = "CTCCCTCTTCTTCTTTTTCATCCTTCTGTCTCTTTGTTTCTGAGCTTTCCTGTCTTTCCTTTTTTCTGAGAGATTCAAAGCCTCCACGACTCTGTTTCCCCCGTCC-CTTCTGAATTTAATTTGCACTAAGTCATTTGCACTGGTTGGAGTTGTGG";
	//read_align = "CTTTCCCTTCGTCTTTTGCTTCTTTTTGT-TTTTTGTTTGTC--CTTTGGTTTCTT-CTTTTTTTTTTTTTTTTTCTGTTCTTCTTTCCCGCTCTCTCCCCTGTTCTCTTCCCTTTCCCTTTTCTTCTTCGTTCGTTGTTGTTGTTGTTGTTGTTG";
	if ((nt >= ntMin && nt <= ntMax))
	{
		string seqtmp = ref_seq_src;
		j = 0;
		for (i = 0; i < nt; i++)
		{
			if (ref_align[i] == '-')
			{
				curchar = '-';
				hashvalue = triplex_score(curchar, read_align[i], Para);
				j++;
			}
			else
			{
				curchar = ref_align_src[j];
				hashvalue = triplex_score(curchar, read_align[i], Para);
//			    cout << curchar << read_align[i] << hashvalue <<endl;
				j++;
			}
			if ((curchar == prechar) && curchar == 'T')
			{
				tri_score = tri_score - prescore + penaltyT;
				hashvalue = penaltyT;
			}
			if ((curchar == prechar) && curchar == 'C')
			{
				tri_score = tri_score - prescore + penaltyC;
				hashvalue = penaltyC;
			}
			prescore = hashvalue;
			if(ref_align[i] != '-')
			prechar = curchar;
			tri_score += hashvalue;
			// TODO: this may be a bug: when hashvalue < 0?
		}
//		cout << tri_score << " is tri_score" << endl;
//		cout << tri_score << " is tri_score" << endl;
//		cout << nt << " is nt" << endl;
		tri_score = tri_score / nt;


	}
//	cout << atr.score << " is sw_score " << atr.stari << " " << atr.endi <<
//	    " is starti endi " << atr.starj << " " << atr.endj <<
//	    " is starj endj " << tri_score << " is tri_score" << endl;
	//cout << nt << " is nt" << endl;

	// tri_score is a little different, maybe due to the different alignment...
	score = (float)alignment.sw_score;
	// Note that when triplex_list is gathered, we need to gather dnaStartPos.
	// Begin to generate triplex_list.
	if ((Para > 0 && strand == 1) || (Para < 0 && strand==0)) {
		refStart = ref_seq.size() - alignment.ref_end - 1;
		refEnd = ref_seq.size() - alignment.ref_begin - 1;
	}
	else {
		refStart = alignment.ref_begin + 1;
		refEnd = alignment.ref_end + 1 ;
	}
	//for reversed DNAseq ,its position is transfered
	//if (Para < 0) reverseSeq(ref_align_src);
	struct triplex fullTriplex;
	fullTriplex = triplex(alignment.query_begin+1, alignment.query_end+1,
		refStart + dnaStartPos, refEnd + dnaStartPos,
		strand, Para, rule, nt, score, identity, tri_score,
		read_align, ref_align_src, 0, 0, 0, 0);
	// filter fullTriplex based on ntMin.
	if (nt >= ntMin)
	{
		triplex_list.push_back(fullTriplex);
	}
	//break; //DEBUG.
	ref_align.clear();
	read_align.clear();
	match = 0;
	mis_match = 0;
	nt = 0;
	score = 0.0;
	identity = 0.0;
	tri_score = 0.0;
	// initialize variables.
}

void getAlignment(StripedSmithWaterman::Alignment alignment,
	string ref_seq,
	string read_seq,
	string ref_seq_src,
	const int8_t* table,
	string &ref_align,
	string &read_align,
	string &ref_align_src)
{
	// Use this function to get alignment of two sequences.
	int cLoop;
	std::vector<uint32_t> tmpCigar;
	for (cLoop = 0; cLoop < alignment.cigar.size(); cLoop++)
	{
		tmpCigar.push_back(alignment.cigar[cLoop]);
	}
	if (tmpCigar.size() > 0)
	{
		// begin to generate alignment.
		int32_t c = 0, left = 0, e = 0, qb = alignment.ref_begin, pb = alignment.query_begin;//need to CHECK
		uint32_t i;
		uint32_t j;
		// seems that we need to determine S???
		//int w = 0, y = 0;
		//for(w = 0; w < alignment.cigar_string.length(); w++)
		//{
		//	if(alignment.cigar_string[w] == 'S')
		//	{
		//		y = w;
		//		break;
		//	}
		//}
		//string ref_align;
		//string read_align;
		while (e < tmpCigar.size() || left > 0)
		{
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			//fprintf(stdout, "Target: %8d		", q + 1);
			// DEBUG.
			//e = y + 1;
			for (c = e; c < tmpCigar.size(); ++c)
			{
				char letter = cigar_int_to_op(tmpCigar[c]);
				uint32_t length = cigar_int_to_len(tmpCigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; ++i)
				{
					if (letter == 'I')
					{
						//fprintf(stdout, "-");
						ref_align = ref_align + "-";
						ref_align_src = ref_align_src + '-';
					}
					else
					{
						//fprintf(stdout, "%c", ref_seq[q]);
						ref_align = ref_align + ref_seq[q];
						ref_align_src = ref_align_src + ref_seq_src[q];
						++q;
					}
					++count;
					if (count == 60) goto step2;
				}
			}//for c = e.

		step2:
			//fprintf(stdout, "		%d\n										", q);
			q = qb;
			count = 0;
			for (c = e; c < tmpCigar.size(); ++c)
			{
				char letter = cigar_int_to_op(tmpCigar[c]);
				uint32_t length = cigar_int_to_len(tmpCigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; ++i)
				{
					if (letter == 'M')
						//if (letter == '=')
					{
						if (table[(int)ref_seq[q]] == table[(int)read_seq[p]])
						{
							//fprintf(stdout, "|");
						}
						else
						{
							//fprintf(stdout, "*");
						}
						++q;
						++p;
					}
					else
					{
						//fprintf(stdout, "*");
						if (letter == 'I')
						{
							++p;
						}
						else
						{
							++q;
						}
					}
					++count;
					if (count == 60)
					{
						qb = q;
						goto step3;
					}
				}
			}// for c = e.
		step3:
			p = pb;
			//fprintf(stdout, "\nQuery:	%8d		", p + 1);
			count = 0;
			for (c = e; c < tmpCigar.size(); ++c)
			{
				char letter = cigar_int_to_op(tmpCigar[c]);
				uint32_t length = cigar_int_to_len(tmpCigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; i++)
				{
					if (letter == 'D')
					{
						//fprintf(stdout, "-");
						read_align = read_align + "-";
					}
					else
					{
						//fprintf(stdout, "%c", read_seq[p]);
						read_align = read_align + read_seq[p];
						++p;
					}
					++count;
					if (count == 60)
					{
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;

					}
				}
			}// for c = e.
			e = c;
			left = 0;
		end:
			//fprintf(stdout, "		%d\n\n", p);
			j = 0;
			//2021-09-16 23:04:55: we don't need to print alignment, we 
			// just need to get alignment sequence.
			//cout << ref_align << " is ref_align" << endl;
			//cout << read_align << " is read_align" << endl;
		}
	}
//	cout << ref_align << "is ref_align" << endl;
//	cout << read_align << "is read_align" << endl;
//	cout << ref_align_src << "is ref_align_src" << endl;
}

