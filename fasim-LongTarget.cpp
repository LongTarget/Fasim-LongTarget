/*
* 2021-09-25 21:38:09: This is a new version of LongTarget which contains some
* new features:
*	 1) now both Sim and fastSIM is available.
*	 2) threshold is determined based on maxScore of the DNA sequence.
* Need to check the performance of fastSim.
* Note that cutSequence has been moved to fastSim.h.
*/
/*
* 2021-12-22 21:38:09: This is a new version of LongTarget which contains some
* new features:
*	 1) TT penalty is operated on sequence with no gaps.
*	 2) threshold is determined based on 0.8*average_maxScore of the 48 DNA transformed sequences.
* Need to check the performance of fastSim.
* Note that cutSequence has been moved to fastSim.h.
*/
#include "fastSim.h"
using namespace std;
struct lgInfo
{
	lgInfo() {};
	lgInfo(const string &s1, const string &s2, const string &s3, const string &s4,
		const string &s5, const string &s6, int s7, const string &s8) :
		lncName(s1), lncSeq(s2), species(s3), dnaChroTag(s4), fileName(s5),
		dnaSeq(s6), startGenome(s7), resultDir(s8) {};
	string lncName;
	string lncSeq;
	string species;
	string dnaChroTag;
	string fileName;
	string dnaSeq;
	int startGenome;
	string resultDir;
};
struct para
{
	string file1path;
	string file2path;
	string outpath;
	int rule;
	int cutLength;
	int strand;
	int overlapLength;
	int minScore;
	bool detailOutput;
	int ntMin;
	int ntMax;
	float scoreMin;
	float minIdentity;
	float minStability;
	int penaltyT;
	int penaltyC;
	int cDistance;
	int cLength;
	bool doFastSim;
	// 2021-09-25 22:21:08: use this parameter to determine do fastSIM or not.
};
struct axis
{
	axis(int n1 = 0, int n2 = 0) :
		triplexnum(n1), neartriplex(n2) {};
	int triplexnum;
	int neartriplex;
};
//void cutSequence(string& seq, vector<string>& seqsVec, vector<int>& seqsStartPos, int& cutLength, int& overlapLength,int& cut_num);
void show_help();
void initEnv(int argc, char * const *argv, struct para &paraList);
void LongTarget(struct para &paraList, string rnaSequence, string dnaSequence,
	vector<struct triplex> &sort_triplex_list);

bool comp(const triplex &a, const triplex &b);
string getStrand(int reverse, int strand);
int same_seq(string &w_str);
void printResult(string &species, struct para paraList, string &lncName,
	string &dnaFile, vector<struct triplex> &sort_triplex_list,
	string &chroTag, string &dnaSequence, int start_genome,
	string &c_tmp_dd, string &c_tmp_length, string &resultDir);

string readDna(string dnaFileName, string &species, string &chroTag,
	string &startGenome);
string readRna(string rnaFileName, string &lncName);
void cluster_triplex(int dd, int length, vector<struct triplex>& triplex_list, map<size_t, size_t> class1[], map<size_t, size_t> class1a[], map<size_t, size_t> class1b[], int class_level);
void print_cluster(int c_level, map<size_t, size_t> class1[], int start_genome, string &chro_info, int dna_size, string &rna_name, int distance, int length, string &outFilePath, string &c_tmp_dd, string &c_tmp_length, vector<struct tmp_class> &w_tmp_class);
int main(int argc, char* const* argv)
{
	struct para paraList;
	vector<struct	lgInfo>	lgList;
	initEnv(argc, argv, paraList);
	char c_dd_tmp[10];
	char c_length_tmp[10];
	int c_loop_tmp = 0;
	string c_tmp_dd;
	string c_tmp_length;
	sprintf(c_dd_tmp, "%d", paraList.cDistance);
	sprintf(c_length_tmp, "%d", paraList.cLength);
	for (c_loop_tmp = 0; c_loop_tmp < strlen(c_dd_tmp); c_loop_tmp++)
	{
		c_tmp_dd += c_dd_tmp[c_loop_tmp];
	}
	for (c_loop_tmp = 0; c_loop_tmp < strlen(c_length_tmp); c_loop_tmp++)
	{
		c_tmp_length += c_length_tmp[c_loop_tmp];
	}
	string lncName;
	string lncSeq;
	string species;
	string dnaChroTag;
	string fileName;
	string dnaSeq;
	string resultDir;
	string startGenomeTmp;
	long startGenome;
	clock_t start, end;
	float cpu_time;
	start = clock();
	dnaSeq = readDna(paraList.file1path, species, dnaChroTag, startGenomeTmp);
	startGenome = atoi(startGenomeTmp.c_str());
	lncSeq = readRna(paraList.file2path, lncName);
	fileName = paraList.file1path.substr(0, paraList.file1path.size() - 3);
	lncName.erase(remove(lncName.begin(), lncName.end(), '\r'), lncName.end());
	lncName.erase(remove(lncName.begin(), lncName.end(), '\n'), lncName.end());
	resultDir = paraList.outpath;
	struct lgInfo algInfo;
	algInfo = lgInfo(lncName, lncSeq, species, dnaChroTag, fileName, dnaSeq,
		startGenome, resultDir);
	lgList.push_back(algInfo);
	int i = 0;
	vector<struct triplex> sort_triplex_list;
	LongTarget(paraList, lgList[i].lncSeq, lgList[i].dnaSeq, sort_triplex_list);
	printResult(lgList[i].species, paraList, lgList[i].lncName,
		lgList[i].fileName, sort_triplex_list, lgList[i].dnaChroTag,
		lgList[i].dnaSeq, lgList[i].startGenome, c_tmp_dd, c_tmp_length,
		lgList[i].resultDir);
	end = clock();
	cout << "finished normally" << endl;
	cpu_time = ((float)(end - start)) / CLOCKS_PER_SEC;
	fprintf(stderr, "CPU time: %f seconds\n", cpu_time);
	return 0;
}

//2021-09-25 21:40:22: cutSequence has been moved to fastSim.h.
//void cutSequence(string& seq, vector<string>& seqsVec, vector<int>& seqsStartPos, int& cutLength, int& overlapLength,int &cut_num)
//{
//	unsigned int pos=0;
//	int tmpa=0;
//	int tmpb=0;
//	seqsVec.clear();seqsStartPos.clear();
//	string cutSeq;
//	while(pos<seq.size())
//	{
//		cutSeq = seq.substr(pos,cutLength);
//		seqsVec.push_back(cutSeq);
//		seqsStartPos.push_back(pos);
//		pos += cutLength;
//		pos -= overlapLength;
//		tmpa++;
//	}
//	cut_num=tmpa;
//}
string readRna(string rnaFileName, string &lncName)
{
	ifstream rnaFile;
	string tmpRNA;
	string tmpStr;
	rnaFile.open(rnaFileName.c_str());
	getline(rnaFile, tmpStr);
	int i = 0;
	string tmpInfo;
	for (i = 0; i < tmpStr.size(); i++)
	{
		if (tmpStr[i] == '>')
		{
			continue;
		}
		tmpInfo = tmpInfo + tmpStr[i];
	}
	lncName = tmpInfo;
	cout << lncName << endl;
	while (getline(rnaFile, tmpStr))
	{
		tmpRNA = tmpRNA + tmpStr;
	}
	return tmpRNA;
}

string readDna(string dnaFileName, string &species, string &chroTag,
	string &startGenome)
{
	ifstream dnaFile;
	string tmpDNA;
	string tmpStr;
	dnaFile.open(dnaFileName.c_str());
	getline(dnaFile, tmpStr);
	int i = 0;
	int j = 0;
	string tmpInfo;
	for (i = 0; i < tmpStr.size(); i++)
	{
		if (tmpStr[i] == '>')
		{
			continue;
		}
		if (tmpStr[i] == '|' && j == 0)
		{
			species = tmpInfo;
			j++;
			tmpInfo.clear();
			continue;
		}
		if (tmpStr[i] == '|' && j == 1)
		{
			chroTag = tmpInfo;
			j++;
			tmpInfo.clear();
			continue;
		}
		if (tmpStr[i] == '-' && j == 2)
		{
			startGenome = tmpInfo;
			tmpInfo.clear();
			continue;
		}
		tmpInfo = tmpInfo + tmpStr[i];
	}
	cout << species << endl;
	cout << chroTag << endl;
	cout << startGenome << endl;
	while (getline(dnaFile, tmpStr))
	{
		tmpDNA = tmpDNA + tmpStr;
	}
	return tmpDNA;
}

void initEnv(int argc, char * const *argv, struct para &paraList)
{
	const char* optstring = "f:s:r:O:c:m:t:i:S:o:y:z:Y:Z:h:D:E:dF";
	struct option long_options[] = {
		{"f1", required_argument, NULL, 'f'},
		{"f2", required_argument, NULL, 's'},
		{"ni", required_argument, NULL, 'y'},
		{"na", required_argument, NULL, 'z'},
		{"pc", required_argument, NULL, 'Y'},
		{"pt", required_argument, NULL, 'Z'},
		{"ds", required_argument, NULL, 'D'},
		{"lg", required_argument, NULL, 'E'},
		{0, 0, 0, 0}
	};
	paraList.file1path = "./";
	paraList.file2path = "./";
	paraList.outpath = "./";
	paraList.rule = 0;
	paraList.cutLength = 5000;
	paraList.strand = 0;
	paraList.overlapLength = 100;
	paraList.minScore = 0;
	paraList.detailOutput = false;
	paraList.ntMin = 20;
	paraList.ntMax = 100000;
	paraList.scoreMin = 0.0;
	paraList.minIdentity = 60.0;
	paraList.minStability = 1;
	paraList.penaltyT = -1000;
	paraList.penaltyC = 0;
	paraList.cDistance = 15;
	paraList.cLength = 50;
	paraList.doFastSim = true;
	int opt;
	if (argc == 1)
	{
		show_help();
	}
	while ((opt = getopt_long_only(argc, argv, optstring, long_options, NULL))
		!= -1)
	{
		switch (opt)
		{
		case 'f':
			paraList.file1path = optarg;
			break;
		case 's':
			paraList.file2path = optarg;
			break;
		case 'r':
			paraList.rule = atoi(optarg);
			break;
		case 'O':
			paraList.outpath = optarg;
			break;
		case 'c':
			paraList.cutLength = atoi(optarg);
			break;
		case 'm':
			paraList.minScore = atoi(optarg);
			break;
		case 't':
			paraList.strand = atoi(optarg);
			break;
		case 'd':
			paraList.detailOutput = true;
			break;
		case 'i':
			paraList.minIdentity = atoi(optarg);
			break;
		case 'S':
			paraList.minStability = atoi(optarg);
			break;
		case 'y':
			paraList.ntMin = atoi(optarg);
			break;
		case 'z':
			paraList.ntMax = atoi(optarg);
			break;
		case 'Y':
			paraList.penaltyC = atoi(optarg);
			break;
		case 'Z':
			paraList.penaltyT = atoi(optarg);
			break;
		case 'o':
			paraList.overlapLength = atoi(optarg);
			break;
		case 'h':
			show_help();
			break;
		case 'D':
			paraList.cDistance = atoi(optarg);
			break;
		case 'E':
			paraList.cLength = atoi(optarg);
			break;
		case 'F':
			paraList.doFastSim = optarg;
		}
	}
}

void LongTarget(struct para &paraList, string rnaSequence, string dnaSequence,
	vector<struct triplex> &sort_triplex_list)
{
	vector< string> dnaSequencesVec;
	vector< int> dnaSequencesStartPos;
	int cut_num = 0;
	cutSequence(dnaSequence, dnaSequencesVec,
		dnaSequencesStartPos, paraList.cutLength,
		paraList.overlapLength, cut_num);
	vector<struct triplex> triplex_list;
	float Identity = 0.0;
	double w_t1, w_t2;
	clock_t time1;
	int minScore = 0,minscore;
	int ji;
	string seqrev;
	string seq2;
	time1 = clock();

	for (int i = 0; i < dnaSequencesVec.size(); i++)
	{
		long dnaStartPos = dnaSequencesStartPos[i];
		cout << "dnaPos = " << dnaStartPos << endl;
		string seq1 = dnaSequencesVec[i];
		if (same_seq(seq1))
		{
			continue;
		}
//		for(int j=1;j<3;j++){
//            if(j==1) {
//                for(int k = 0; k<6;k++){
//                     seq2 = transferString(seq1, 0, 1, k + 1);
//                     minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule);
//                     cout<<"____________________"<<minscore<<endl;
//                     minScore += minscore;
//                     seq2 = transferString(seq1, 1, 1, k + 1);
//                    //cout<<"____________________"<<j<<endl;
//                     minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule);
//                                          cout<<"____________________"<<minscore<<endl;
//                     minScore += minscore;
//                }
//            }
//            else{
//                 for(int k = 0; k<18;k++){
//                     seq2 = transferString(seq1, 0, -1, k + 1);
//                     //cout<<"____________________"<<j<<endl;
//                     minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule);
//                                          cout<<"____________________"<<minscore<<endl;
//                     minScore += minscore;
//                     seq2 = transferString(seq1, 1, -1, k + 1);
//                    //cout<<"____________________"<<j<<endl;
//                     minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule);
//                                          cout<<"____________________"<<minscore<<endl;
//                     minScore += minscore;
//                }
//            }
//        }
//        minScore /= 48.0;
//        minScore *= 0.8;
//        cout<<"minScore____________________"<<minScore<<endl;
		if (paraList.strand >= 0)
		{
			if (paraList.rule == 0)
			{
				for (int j = 0; j < 6; j++)
				{
					string seq2 = transferString(seq1, 0, 1, j + 1);
					//cout<<"____________________"<<j<<endl;
					if (paraList.doFastSim)
					{
					    //minScore *= 0.8;
					    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
					    //minScore = minScore>minscore?minScore:minscore;
					    minScore = minscore;
						fastSIM(rnaSequence, seq2, seq1, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 0, 1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}
					else
					{
					    //minScore *= 0.8;
			            minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
			            minScore = minscore;
						SIM(rnaSequence, seq2, seq1, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 0, 1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}

					seq2 = transferString(seq1, 1, 1, j + 1);
					reverseSeq(seq2);
					seqrev = seq1;
					complement(seqrev);
					reverseSeq(seqrev);
					if (paraList.doFastSim)
					{
				        minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
                        //minScore = minScore>minscore?minScore:minscore;
                        minScore = minscore;
						fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 1, 1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}
					else
					{
					    //minScore *= 0.8;
				        minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				        minScore = minscore;
						SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 1, 1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}
				}
			}
			if (paraList.rule > 0 && paraList.rule < 7)
			{
				string seq2 = transferString(seq1, 0, 1, paraList.rule);

				if (paraList.doFastSim)
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
                    //minScore = minScore>minscore?minScore:minscore;
                    minScore = minscore;
					fastSIM(rnaSequence, seq2, seq1, dnaStartPos, minScore, 5, -4,
						-12, -4, triplex_list, 0, 1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}
				else
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				    minScore = minscore;
					SIM(rnaSequence, seq2, seq1, dnaStartPos, minScore, 5, -4,
						-12, -4, triplex_list, 0, 1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}

				seq2 = transferString(seq1, 1, 1, paraList.rule);
				reverseSeq(seq2);
				seqrev = seq1;
				complement(seqrev);
				reverseSeq(seqrev);
				if (paraList.doFastSim)
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
                    //minScore = minScore>minscore?minScore:minscore;
                    minScore = minscore;
					fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 1, 1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}
				else
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				    minScore = minscore;
					SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 1, 1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}
			}
		}
		if (paraList.strand <= 0)
		{
			if (paraList.rule == 0)
			{
				for (int j = 0; j < 18; j++)
				{
					string seq2 = transferString(seq1, 1, -1, j + 1);
					seqrev = seq1;
					complement(seqrev);
					//cout << "______1____" << j+1 << "__________" << minScore <<endl;
					if (paraList.doFastSim)
					{
					    //minScore *= 0.8;
				        minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
					    //minScore = minScore>minscore?minScore:minscore;
					    minScore = minscore;
						fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 1, -1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}
					else
					{
					    //minScore *= 0.8;
					    //cout<<j+1<<"______minScore____________________"<<minScore<<endl;
					    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
					    minScore = minscore;
						SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 1, -1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}

					seq2 = transferString(seq1, 0, -1, j + 1);
					reverseSeq(seq2);
					seqrev = seq1;
					reverseSeq(seqrev);
					//cout << "_____0_____" << j+1 << "__________" << minScore <<endl;
					//cout << seq2 <<endl;
					if (paraList.doFastSim)
					{
					    //minScore *= 0.8;
				        minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
					    //minScore = minScore>minscore?minScore:minscore;
					    minScore = minscore;
						fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 0, -1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}
					else
					{
					    //minScore *= 0.8;
					    //cout<<j+1<<"______minScore____________________"<<minScore<<endl;
					    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
					    minScore = minscore;
						SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4,
							-12, -4, triplex_list, 0, -1, j + 1, paraList.ntMin,
							paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
					}

				}
			}
			else
			{
				string seq2 = transferString(seq1, 1, -1, paraList.rule);
				seqrev = seq1;
				complement(seqrev);
				if (paraList.doFastSim)
				{
				   // minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				   	//minScore = minScore>minscore?minScore:minscore;
				   	minScore = minscore;
					fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 1, -1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}
				else
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				    minScore = minscore;
					SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 1, -1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}

				seq2 = transferString(seq1, 0, -1, paraList.rule);
				reverseSeq(seq2);
				//minScore = calc_score_once(rnaSequence, seq2, dnaStartPos,
					//paraList.rule);
				seqrev = seq1;
				reverseSeq(seqrev);
				if (paraList.doFastSim)
				{
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				   	//minScore = minScore>minscore?minScore:minscore;
				   	minScore = minscore;
					fastSIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 0, -1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}
				else
				{
				    //minScore *= 0.8;
				    minscore = calc_score_once(rnaSequence, seq2, dnaStartPos,paraList.rule)*0.8;
				    minScore = minscore;
					SIM(rnaSequence, seq2, seqrev, dnaStartPos, minScore, 5, -4, -12,
						-4, triplex_list, 0, -1, paraList.rule, paraList.ntMin,
						paraList.ntMax, paraList.penaltyT, paraList.penaltyC);
				}

			}
		}
	}
	clock_t time2;
	time2 = clock();
	for (int i = 0; i < triplex_list.size(); i++)
	{
		triplex atr = triplex_list[i];
		//cout << atr.stari <<'\t'<<atr.endi <<'\t'<< atr.starj <<'\t'<<atr.endj << atr.strand <<'\t'<<atr.reverse<<'\t'<< atr.score <<'\t'<<atr.nt<<'\t'<<atr.identity<<'\t'<<atr.nt<<'\t'<<atr.tri_score<<endl;

		if (atr.score >= paraList.scoreMin && atr.identity >= paraList.minIdentity
			&& atr.tri_score >= paraList.minStability && atr.nt >= paraList.cLength)
		{
			sort_triplex_list.push_back(atr);
		}
	}
}

void cluster_triplex(int dd, int length, vector<struct triplex>& triplex_list, map<size_t, size_t> class1[], map<size_t, size_t> class1a[], map<size_t, size_t> class1b[], int class_level)
{
	int i, j;
	int find = 0;
	map<size_t, struct axis> axis_map;
	int max_neartriplexnum = 0, max_pos = 0;
	int middle = 0;
	int count = 0;
	for (vector<struct triplex>::iterator it = triplex_list.begin(); it != triplex_list.end(); it++)
	{
	    //cout<< it->stari <<"-------"<<it->endi<<"  "<< it->nt <<endl;
		if (it->nt > length)
		{
			count++;
			middle = (int)((it->stari + it->endi) / 2);
			it->middle = middle;
			it->motif = 0;
			axis_map[middle].triplexnum++;

			for (i = -dd; i <= dd; i++)
			{
				if (i > 0)
				{
					axis_map[middle + i].neartriplex = axis_map[middle + i].neartriplex + (dd - i);
				}
				else if (i < 0)
				{
					axis_map[middle + i].neartriplex = axis_map[middle + i].neartriplex + (dd + i);
				}
				else
				{
				}
				if (axis_map[middle].triplexnum > 0)
				{
				    //cout<< middle+i << " hit>1  " << axis_map[middle + i].neartriplex <<endl;
					if (axis_map[middle + i].neartriplex > max_neartriplexnum)
					{
						max_neartriplexnum = axis_map[middle + i].neartriplex;
						max_pos = middle + i;
						find = 1;
					}
				}
			}
			it->neartriplex = axis_map[middle].neartriplex;
		}
	}
	int theclass = 1;
	while (find)
	{
		for (i = max_pos - dd; i <= max_pos + dd; i++)
		{
			for (vector<struct triplex>::iterator it = triplex_list.begin(); it != triplex_list.end(); it++)
			{
				if (it->middle == i && it->motif == 0)
				{
					it->motif = theclass;
					it->center = max_pos;
					if (theclass > class_level)
					{
						continue;
					}
					if (it->endj > it->starj)
						for (j = it->starj; j < it->endj; j++)
						{
							class1[theclass][j]++;
							class1a[theclass][j]++;
						}
					else
						for (j = it->endj; j < it->starj; j++)
						{
							class1[theclass][j]++;
							class1b[theclass][j]--;
						}
				}
			}
			//cout<<"axis_map.erase  "<< axis_map[i].neartriplex << "  pos "<< i <<endl;
			axis_map.erase(i);
		}
		max_neartriplexnum = 0;
		find = 0;
		/*for (map<size_t, struct axis>::iterator it = axis_map.begin(); it != axis_map.end(); it++)
		{
			if (it->second.neartriplex >= max_neartriplexnum && it->second.triplexnum > 0)
			{
			    cout<< it->second.neartriplex <<"   " << max_neartriplexnum <<"   max_neartriplexnum  "<<max_pos<<"  max_pos  "<<endl;
				max_neartriplexnum = it->second.neartriplex;
				max_pos = it->first;
				find = 1;
			}
		}*/
		for (i = 0 ; i<axis_map.size(); i++)
		{
		    //if(axis_map[i].neartriplex>0)
		    //cout<< "  left_position  "<<axis_map[i].neartriplex  << "     pos    " << i <<endl;
			if (axis_map[i].neartriplex > max_neartriplexnum)
			{
			    //cout<< axis_map[i].neartriplex <<"   " << max_neartriplexnum <<"   max_neartriplexnum  "<<max_pos<<"  max_pos  "<<endl;
				max_neartriplexnum = axis_map[i].neartriplex;
				max_pos = i;
				find = 1;
			}
		}
		++theclass;
	}
}


void print_cluster(int c_level, map<size_t, size_t> class1[], int start_genome, string &chro_info, int dna_size, string &rna_name, int distance, int length, string &outFilePath, string &c_tmp_dd, string &c_tmp_length, vector<struct tmp_class> &w_tmp_class)
{
	struct tmp_class a_tmp_class;
	char c_level_tmp[3];
	cout << c_level_tmp << c_level << endl;
	sprintf(c_level_tmp, "%d", c_level);
	string c_tmp_level;
	int c_level_loop = 0;
	for (c_level_loop = 0; c_level_loop < strlen(c_level_tmp); c_level_loop++)
	{
		c_tmp_level += c_level_tmp[c_level_loop];
	}
	string class_name = outFilePath.substr(0, outFilePath.size() - 17) + "-TFOclass" + c_tmp_level;
	ofstream outfile(class_name.c_str(), ios::trunc);
	int map_tmp0 = 0, map_tmp1 = 0, map_tmp2 = 0, map_tmp3 = 0, map_count = 0, map_count1 = 0;
	int map_first1 = 0, map_second1 = 0;
	int map_first0 = 0, map_second0 = 0;
	int if_map1 = 0, if_map2 = 0, if_map3 = 0, if_map4 = 0;
	int if_map_flag = 0;
	outfile << "browser position " << chro_info << ":" << start_genome << "-" << start_genome + dna_size << endl;
	outfile << "browser hide all" << endl;
	outfile << "browser pack refGene encodeRegions" << endl;
	outfile << "browser full altGraph" << endl;
	outfile << "# 300 base wide bar graph, ausoScale is on by default == graphing" << endl;
	outfile << "# limits will dynamically change to always show full range of data" << endl;
	outfile << "# in viewing window, priority = 20 position this as the second graph" << endl;
	outfile << "# Note, zero-relative, half-open coordinate system in use for bedGraph format" << endl;
	outfile << "track type=bedGraph name='" << rna_name << " TTS (" << c_level << ")' description='" << distance << "-" << length << "' visibility=full color=200,100,0 altColor=0,100,200 priority=20" << endl;
	int final_genome = 0;
	for (map<size_t, size_t>::iterator it = class1[c_level].begin(); it != class1[c_level].end(); it++)
	{
		final_genome = it->first + start_genome;
	}
	for (map<size_t, size_t>::iterator it = class1[c_level].begin(); it != class1[c_level].end(); )
	{
		map_first0 = it->first;
		map_tmp1 = it->first;
		map_tmp2 = it->second;
		if ((it->first + start_genome) == final_genome || it == class1[c_level].end())
		{
			a_tmp_class = tmp_class(map_first0 + start_genome - 1, map_tmp1 + start_genome, map_tmp2, 0, 0);
			w_tmp_class.push_back(a_tmp_class);
			break;
		}
		it++;
		while (abs((long)(it->first - map_tmp1)) == 1 && (it->second == map_tmp2))
		{
			if ((it->first + start_genome) == final_genome)
			{
				break;
			}
			map_tmp1 = it->first;
			map_tmp2 = it->second;
			it++;
		}
		if (map_count == 0)
		{
			a_tmp_class = tmp_class(map_first0 + start_genome - 2, map_tmp1 + start_genome, map_tmp2, 0, 0);
			w_tmp_class.push_back(a_tmp_class);
			map_count++;
		}
		else
		{
			a_tmp_class = tmp_class(map_first0 + start_genome - 1, map_tmp1 + start_genome, map_tmp2, 0, 0);
			w_tmp_class.push_back(a_tmp_class);
		}
		if (abs((long)(it->first - map_tmp1)) != 1)
		{
			a_tmp_class = tmp_class(map_tmp1 + start_genome, it->first + start_genome - 1, 0, 0, 0);
			w_tmp_class.push_back(a_tmp_class);

		}
	}
	int w_class_loop = 0;
	for (w_class_loop = 0; w_class_loop < w_tmp_class.size(); w_class_loop++)
	{
		tmp_class btc = w_tmp_class[w_class_loop];
		outfile << chro_info << "\t" << btc.genome_start << "\t" << btc.genome_end << "\t" << btc.signal_level << endl;

		/*tmp_class ctc = w_tmp_class[w_class_loop + 1];
		if (btc.genome_start == final_genome)
		{
			break;
		}
		if (w_class_loop + 1 == w_tmp_class.size())
		{
		}
		if (btc.genome_start == ctc.genome_start)
		{
			if (1)
			{
				outfile << chro_info << "\t" << btc.genome_start << "\t" << ctc.genome_end << "\t" << ctc.signal_level << endl;
			}
			w_class_loop += 1;
		}
		else
		{
			outfile << chro_info << "\t" << btc.genome_start << "\t" << btc.genome_end << "\t" << btc.signal_level << endl;

		}*/
	}
}

void printResult(string &species, struct para paraList, string &lncName, string &dnaFile, vector<struct triplex> &sort_triplex_list, string &chroTag, string &dnaSequence, int start_genome, string &c_tmp_dd, string &c_tmp_length, string &resultDir)
{
	vector<struct tmp_class> w_tmp_class;
	string pre_file2 = resultDir + "/" + species + "-" + lncName;
	string outFilePath = pre_file2 + "-fastSim-TFOsorted";
	ofstream outFile(outFilePath.c_str(), ios::trunc);
	outFile << "QueryStart\t" << "QueryEnd\t" << "StartInSeq\t" << "EndInSeq\t" << "Direction\t" << "StartInGenome\t" << "EndInGenome\t" << "MeanStability\t" << "MeanIdentity(%)\t" << "Strand\t" << "Rule\t" << "Score\t" << "Nt(bp)\t" << "Class\t" << "MidPoint\t" << "Center\t" << "TFO sequence" << "TTS sequence"<< endl;
	map<size_t, size_t> class1[6], class1a[6], class1b[6];
	int class_level = 5;
	cluster_triplex(paraList.cDistance, paraList.cLength, sort_triplex_list, class1, class1a, class1b, class_level);
	sort(sort_triplex_list.begin(), sort_triplex_list.end(), comp);
	for (int i = 0; i < sort_triplex_list.size(); i++)
	{
		triplex atr = sort_triplex_list[i];
		if (sort_triplex_list[i].motif == 0)
		{
			continue;
		}

//		atr.stri_align.erase(remove(atr.stri_align.begin(), atr.stri_align.end(), '-'), atr.stri_align.end());
		if (atr.starj < atr.endj)
			outFile << atr.stari << "\t" << atr.endi << "\t" << atr.starj << "\t" << atr.endj << "\t" << "R\t" << atr.starj + start_genome - 1 << "\t" << atr.endj + start_genome - 1 << "\t" << atr.tri_score << "\t" << atr.identity << "\t" << getStrand(atr.reverse, atr.strand) << "\t" << atr.rule << "\t" << atr.score << "\t" << atr.nt << "\t" << atr.motif << "\t" << atr.middle << "\t" << atr.center << "\t" << atr.stri_align << "\t" << atr.strj_align<< endl;
		else
			outFile << atr.stari << "\t" << atr.endi << "\t" << atr.starj << "\t" << atr.endj << "\t" << "L\t" << atr.endj + start_genome - 1 << "\t" << atr.starj + start_genome - 1 << "\t" << atr.tri_score << "\t" << atr.identity << "\t" << getStrand(atr.reverse, atr.strand) << "\t" << atr.rule << "\t" << atr.score << "\t" << atr.nt << "\t" << atr.motif << "\t" << atr.middle << "\t" << atr.center << "\t" << atr.stri_align << "\t" << atr.strj_align<< endl;
	}
	outFile.close();
	int pr_loop = 0;
	for (pr_loop = 1; pr_loop < 3; pr_loop++)
	{
		print_cluster(pr_loop, class1, start_genome - 1, chroTag, dnaSequence.size(), lncName, paraList.cDistance, paraList.cLength, outFilePath, c_tmp_dd, c_tmp_length, w_tmp_class);
	}
	vector<struct tmp_class>tmpClass;
	tmpClass.swap(w_tmp_class);
	for (pr_loop = 0; pr_loop < 6; pr_loop++)
	{
		class1[pr_loop].clear();
		class1a[pr_loop].clear();
		class1b[pr_loop].clear();
	}
}

bool comp(const triplex &a, const triplex &b)
{
	return a.motif < b.motif;
}
string getStrand(int reverse, int strand)
{
	string Strand;
	if (reverse == 1 && strand == 0)
	{
		Strand = "ParaPlus";
	}
	else if (reverse == 1 && strand == 1)
	{
		Strand = "ParaMinus";
	}
	else if (reverse == -1 && strand == 1)
	{
		Strand = "AntiMinus";
	}
	else if (reverse == -1 && strand == 0)
	{
		Strand = "AntiPlus";
	}
	return Strand;
}

int same_seq(string &w_str)
{
	string A = w_str;
	int i = 0;
	int a = 0, c = 0, g = 0, t = 0, u = 0, n = 0;
	for (i = 0; i < A.size(); i++)
	{
		switch (A[i])
		{
		case 'A':
			a++;
			break;
		case 'C':
			c++;
			break;
		case 'G':
			g++;
			break;
		case 'T':
			t++;
			break;
		case 'U':
			u++;
			break;
		case 'N':
			n++;
			break;
		default:
			cout << "unknown letter" << endl;
			break;
		}
	}
	if (a == A.size())
	{
		return 1;
	}
	else if (c == A.size())
	{
		return 1;
	}
	else if (g == A.size())
	{
		return 1;
	}
	else if (t == A.size())
	{
		return 1;
	}
	else if (u == A.size())
	{
		return 1;
	}
	else if (n == A.size())
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void show_help()
{
	cout << "This is the help page." << endl;
	cout << "options	 Parameters			functions" << endl;
	cout << "f1	 DNA sequence file	used to get the DNA sequence" << endl;
	cout << "f2	 RNA sequence file	used to get the RNA sequence" << endl;
	cout << "r		rules							rules used to construct triplexes.int type.0 is all." << endl;
	cout << "O		Output path				if you define this,output result will be in the path.default is pwd" << endl;
	cout << "c		Cutlength					Cut sequence's length." << endl;
	cout << "m		min_score					Min_score...this option maybe useless.keep it for now." << endl;
	cout << "d		detailoutut				if you choose -d option,it will generate a triplex.detail file which describes the sequence-alignment." << endl;
	cout << "i		identity					 a condition used to pick up triplexes.default is 60.this should be int type such as 60,not 0.6.default is 60." << endl;
	cout << "S		stability					a condition like identity,should be float type such as 1.0.default is 1.0." << endl;
	cout << "ni	 ntmin							triplexes' min length.default is 20." << endl;
	cout << "na	 ntmax							triplexes' max length.default is 100." << endl;
	cout << "pc	 penaltyC					 penalty about GG.default is 0." << endl;
	cout << "pt	 penaltyT					 penalty about AA.default is -1000." << endl;
	cout << "ds	 c_dd							 distance used by cluster function.default is 15." << endl;
	cout << "lg	 c_length					 triplexes' length threshold used in cluster function.default is 50." << endl;
	cout << "F     doFastSim     if true, fastSIM function will be used instead of SIM function." << endl;
	cout << "all parameters are listed.If you want to run a simple example,type ./LongTarget -f1 DNAseq.fa -f2 RNAseq.fa -r 0 will be OK" << endl;
	cout << "any problems or bugs found please send email to us:zhuhao@smu.edu.cn." << endl;
	exit(1);
}
