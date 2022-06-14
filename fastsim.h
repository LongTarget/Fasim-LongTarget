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

struct para
{
	string file1path;
	string file2path;
	string outpath;
	int corenum;
	int rule;
	int cutLength;
	int strand;
	int overlapLength;
	int minScore;
	bool detailOutput;
	bool doFastSim;
	int ntMin;
	int ntMax;
	float scoreMin;
	float minIdentity;
	float minStability;
	int penaltyT;
	int penaltyC;
	int cDistance;
	int cLength;
	// 2021-09-25 22:21:08: use this parameter to determine do fastSIM or not.
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
	int penaltyC, struct para paraList)
{
	int32_t maskLen = 15;
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment alignment;
	struct mycutregion tmpregion;
	vector<struct mycutregion> mycut_list;
	vector<struct triplex> myTriplexList;
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

	for (i = 0; i < finalScoreInfo.size(); i++)
	{
		float Iden = 0.6;
		int cutlength,bestcutregion;
		int myflag=0;
		StripedSmithWaterman::Alignment bestalignment;
		bestalignment.sw_score = 0;
		while(Iden<=1){
            cutlength = (int)(finalScoreInfo[i].score+24)/(9*Iden-4)+1;
            cutlength = finalScoreInfo[i].position - cutlength+1  > 0 ? cutlength:finalScoreInfo[i].position+1;
            smallSeq = strB.substr(finalScoreInfo[i].position - cutlength + 1, cutlength);
//            cout<<finalScoreInfo[i].score<<"\tfinalscore\t"<<finalScoreInfo[i].position<<"\tfinalpos\t"<<"cutlength\t"<<cutlength<<endl;
//            cout<<rule<<"\t"<<strand<<"\t"<<Para<<"\n";
//            cout<<smallSeq<<endl;
            aligner.Align(strA.c_str(), smallSeq.c_str(), smallSeq.size(), filter,&alignment, maskLen);
//          cout<<alignment.sw_score<<"\t"<<alignment.ref_end<<endl;
            if(alignment.sw_score>=finalScoreInfo[i].score) {
                myflag = 1;
                break;
                }//||(alignment.ref_end == cutlength)
            if(alignment.sw_score>bestalignment.sw_score&&alignment.ref_end==cutlength-1){
                bestalignment.sw_score = alignment.sw_score;
                bestalignment.sw_score_next_best = alignment.sw_score_next_best;
                bestalignment.ref_begin = alignment.ref_begin;
                bestalignment.ref_end = alignment.ref_end;
                bestalignment.query_begin = alignment.query_begin;
                bestalignment.query_end = alignment.query_end;
                bestalignment.ref_end_next_best = alignment.ref_end_next_best;
                bestalignment.mismatches = alignment.mismatches;
                bestalignment.cigar_string = alignment.cigar_string;
                bestalignment.cigar = alignment.cigar;
                bestcutregion = cutlength;
                myflag=2;
            }
            Iden += 0.1;
		}
		if(myflag==2){
            alignment.sw_score = bestalignment.sw_score;
			alignment.sw_score_next_best = bestalignment.sw_score_next_best;
			alignment.ref_begin = bestalignment.ref_begin;
			alignment.ref_end = bestalignment.ref_end;
			alignment.query_begin = bestalignment.query_begin;
			alignment.query_end = bestalignment.query_end;
			alignment.ref_end_next_best = bestalignment.ref_end_next_best;
			alignment.mismatches = bestalignment.mismatches;
			alignment.cigar_string = bestalignment.cigar_string;
			alignment.cigar = bestalignment.cigar;
            cutlength = bestcutregion;
        }
        bestalignment.Clear();
//        cout<<alignment.sw_score<<"\t"<<alignment.ref_end<<endl;
        if(alignment.sw_score!=0){
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
		    mycut_list.clear();
        }
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
	    triplex atr = myTriplexList[i];
        if(atr.identity >= paraList.minIdentity && atr.tri_score >= paraList.minStability && atr.nt >= ntMin)
		triplex_list.push_back(atr);
	}
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
	struct triplex fullTriplex;
	fullTriplex = triplex(alignment.query_begin+1, alignment.query_end+1,
		refStart + dnaStartPos, refEnd + dnaStartPos,
		strand, Para, rule, nt, score, identity, tri_score,
		read_align, ref_align_src, 0, 0, 0, 0,0,0,"");
	if (nt >= ntMin)
	{
		triplex_list.push_back(fullTriplex);
	}
	ref_align.clear();
	read_align.clear();
	match = 0;
	mis_match = 0;
	nt = 0;
	score = 0.0;
	identity = 0.0;
	tri_score = 0.0;
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
}

