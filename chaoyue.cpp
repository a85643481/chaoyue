// YinYang_k_means_ori.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include<iostream>
#include<limits>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include<string>
#include<string.h>
#include<math.h>
#include <time.h>
#include<vector>
#include <iterator>
#include <algorithm>
#define DBL_MAX 1.7976931348623158e+308
#define max 100
using namespace std;
void datasgetKmeansYinYang(int colum1, int colum2, int num, int k, int &count, double **point, int flag, double**datas, double***pp, double &count1, int tan, double &localCount, double &noLocalCount, int *listclusterIndex);
void getMeanTwo(int *ind, double **datas, double **arr, double**point, int num, int colum, int k, int *size);
double calObjFuncForGivenCents(double **data, double **V, int data_num, int data_dim, int clu_num);
inline double getDistance(double *a, double *b, int colum);
void readText(string file, int colum, int num, double**data);
void readPoint(string file, int colum, int times, int **points);
void datasgetKmeans(int colum1, int colum2, int num, int k, int &count, int flag, double**datas1, double***pp, double &count1, int tan, double &localCount, double &noLocalCount);
inline void getMean(double**datas, int *data1, double data2[], int size, int colum);
inline bool equalsArr(double *vec1, double *vec2, int colum);
void getArr(int a, int b, double **c);
double min(double *arr, int length);
void bubble_sort(int a[], int n);
int num = 0;
int times = 0;
int dim = 0;
int K = 0;
int x = 10;
string to_String(int n)
{
	int m = n;
	char s[max];
	char ss[max];
	int i = 0, j = 0;
	if (n < 0)// 处理负数
	{
		m = 0 - m;
		j = 1;
		ss[0] = '-';
	}
	while (m>0)
	{
		s[i++] = m % 10 + '0';
		m /= 10;
	}
	s[i] = '\0';
	i = i - 1;
	while (i >= 0)
	{
		ss[j++] = s[i--];
	}
	ss[j] = '\0';
	return ss;
}
int _tmain(int argc, _TCHAR* argv[]){
	string dataFile = "";
	cout << "请输入数据文件名字" << endl;
	cin >> dataFile;
	cout << "请输入数据个数" << endl;
	cin >> num;
	cout << "请输入数据维度" << endl;
	cin >> dim;
	cout << "请输入聚类个数k" << endl;
	cin >> K;
	cout << "请输入运行次数" << endl;
	cin >> times;
	cout << "请输入群组个数k/x中的x值" << endl;
	cin >> x;
	/*x =10;*/
	double **datas = new double*[num];
	*datas = new double[num*dim];
	for (int i = 0; i < num; i++){
		datas[i] = *datas + i*dim;
	}
	int **points = new int*[times];
	for (int i = 0; i < times; i++){
		points[i] = new int[K];
	}

	string str = "D:/Fast KMeans Tested Data Sets";
	time_t t = time(0);
	char ch[64];
	strftime(ch, sizeof(ch), "%Y-%m-%d %H-%M-%S", localtime(&t));
	string path = dataFile + "_c=" + to_String(K) + "_dim=" + to_String(dim) + "_dataNum=" + to_String(num) + "_maxRunTimes=" + to_String(times) + "_yinyang" + "(" + "t=" + to_String((K - 1) / x + 1) + "=" + to_String(K) + "：" + to_String(x) + ")_" + ch + ".txt";
	cout << path << endl;
	dataFile = str + "/" + dataFile + ".txt";
	cout << dataFile << endl;
	readText(dataFile, dim, num, datas);
	string pointFile = "N=" + to_String(num) + "_c=" + to_String(K) + "_Runs=50.txt";
	cout << pointFile << endl;
	pointFile = "D:/initIDXForCenters/" + pointFile;
	readPoint(pointFile, K, times, points);
	ofstream fout(path);
	fout << setw(35) << "runTimes";
	fout << setw(35) << "cluNum";
	fout << setw(35) << "dataNum";
	fout << setw(35) << "dataDim";
	fout << setw(35) << "iterTimes";
	fout << setw(35) << "objFunc";
	fout << setw(35) << "distCallsAll";
	fout << setw(35) << "distCallPerIterPerData";
	fout << setw(35) << "      [减少距离计算次数比例 (所有)]";
	fout << setw(35) << "      [减少距离计算次数比例 (非local)]";
	fout << setw(35) << "             timeAll(毫秒)";
	fout << setw(35) << "                   timePerIter(毫秒)" << endl;;
	double sum = 0;
	double sumLocalCountPer = 0.0;
	double sumNoLocalCountPer = 0.0;
	double sumLocalCount = 0;
	double sumNoLocalCount = 0;
	double sumLocalCountIter = 0;
	double sumNoLocalCountIter = 0;
	double countsum = 0;
	double objsum = 0;
	for (int i = 0; i < times; i++){
		double **point = new double*[K];
		for (int j = 0; j < K; j++){
			point[j] = datas[points[i][j]];
		}
		int count = 0;
		double count1 = 0;

		double begintime, endtime;
		begintime = clock();
		double **pp = NULL;
		double countLocal = 0;
		double countNoLocal = 0;
		datasgetKmeans(K, dim, num, K, count, 0, datas, &pp, count1, x, countLocal, countNoLocal);
		endtime = clock();
		sum = sum + endtime - begintime;
		countsum += count;
		double obj = calObjFuncForGivenCents(datas, pp, num, dim, K);
		delete[] * pp;
		objsum += obj;
		countLocal += count1;
		countNoLocal += count1;
		sumLocalCount += countLocal;
		sumNoLocalCount += countNoLocal;
		sumLocalCountPer += (countLocal / ((double)count*(double)num*(double)K*1.0));
		sumNoLocalCountPer += (countNoLocal / ((double)count*(double)num*(double)K*1.0));
		sumLocalCountIter += (countLocal / count) / num;
		sumNoLocalCountIter += (countNoLocal / count) / num;
		fout << setw(35) << setprecision(20) << i + 1;
		fout << setw(35) << setprecision(20) << K;
		fout << setw(35) << setprecision(20) << num;
		fout << setw(35) << setprecision(20) << dim;
		fout << setw(35) << setprecision(20) << count;
		fout << setw(35) << setprecision(20) << obj;
		fout << setw(35) << setprecision(20) << countLocal;
		fout << setw(35) << setprecision(20) << (countLocal / count) / num;
		fout << setw(35) << setprecision(20) << 1 - countLocal / ((double)count*(double)num*(double)K*1.0);
		fout << setw(35) << setprecision(20) << 1 - countNoLocal / ((double)count*(double)num*(double)K*1.0);
		fout << setw(35) << setprecision(20) << endtime - begintime;
		fout << setw(35) << setprecision(20) << (endtime - begintime) / count << endl;

	}
	fout << "/*********************************************avg********************************************************/" << endl;
	fout << setw(35) << setprecision(20) << "-";
	fout << setw(35) << setprecision(20) << "-";
	fout << setw(35) << setprecision(20) << "-";
	fout << setw(35) << setprecision(20) << "-";
	fout << setw(35) << setprecision(20) << countsum / times;
	fout << setw(35) << setprecision(20) << objsum / times;
	fout << setw(35) << setprecision(20) << sumLocalCount / times;
	fout << setw(35) << setprecision(20) << sumLocalCountIter / times;
	fout << setw(35) << setprecision(20) << 1 - sumLocalCountPer / times;
	fout << setw(35) << setprecision(20) << 1 - sumNoLocalCountPer / times;
	fout << setw(35) << setprecision(20) << sum / times;
	fout << setw(35) << setprecision(20) << sum / countsum;
	delete[] * datas;
	delete[] datas;
	fout.close();
	cout << "程序结束";
	getchar();
	getchar();
	return 0;
}
void readText(string file, int colum, int num, double**data){
	ifstream infile;
	infile.open(file.data());
	assert(infile.is_open());

	for (int i = 0; i < num; i++)
	{

		for (int j = 0; j < colum; j++)
		{
			infile >> data[i][j];
		}
	}

	infile.close();
}
void readPoint(string file, int colum, int times, int **points){
	ifstream infile;
	infile.open(file.data());
	assert(infile.is_open());
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < colum; j++)
		{
			infile >> points[i][j];
		}
	}
	infile.close();
}
inline double getDistance(double *a, double *b, int colum){
	double sum = 0;
	int i;
	for (i = 0; i < colum; i++){
		sum += (a[i] - b[i])*(a[i] - b[i]);
	}
	return sqrt(sum);
}
inline double getDistance1(double *a, double *b, int colum){
	double sum = 0;
	int i;
	for (i = 0; i < colum; i++){
		sum += (a[i] - b[i])*(a[i] - b[i]);
	}
	return sum;
}
void bubble_sort(int a[], int n)
{
	int i, j, temp;
	for (j = 0; j < n - 1; j++)
	for (i = 0; i < n - 1 - j; i++)
	if (a[i] > a[i + 1])
	{
		temp = a[i]; a[i] = a[i + 1]; a[i + 1] = temp;
	}
}
double** getArr(int a, int b){
	double **arr = new double*[a];
	*arr = new double[a*b];
	for (int i = 0; i < a; i++)
	{
		arr[i] = *arr + i*b;
	}
	return arr;
}

int** getIntegerArr(int a, int b){
	int **arr = new int*[a];
	*arr = new int[a*b];
	for (int i = 0; i < a; i++)
	{
		arr[i] = *arr + i*b;
	}
	return arr;
}
//给数组初始化0
void setZeroDouble(double*a, int n){
	for (int i = 0; i < n; i++){
		a[i] = 0;
	}
}
//给数组初始化0
void setZeroInteger(int*a, int n){
	for (int i = 0; i < n; i++){
		a[i] = 0;
	}
}
//拷贝一个数组
void copyDouble(double*a, double*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//拷贝一个数组
void copyInteger(int*a, int*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//拷贝一个数组
void copyBool(bool*a, bool*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//向量除以一个数赋值给另一个向量
void arrDivide(double *a, double b, int length, double *arr){
	for (int i = 0; i < length; i++){
		arr[i] = a[i] / b;
	}
}
void arrDivide(double *a, double b, int length){
	for (int i = 0; i < length; i++){
		a[i] = a[i] / b;
	}
}
//两个向量相加
void arrAdd(double *a, double *b, int length){
	for (int i = 0; i < length; i++){
		a[i] += b[i];
	}
}
void getRamdArr(double ***points, double **datas, int Groupcount, int colum2, int num, int *groupDateSize, int *indexClis, int size){
	for (int j = 0; j < size; j++){
		groupDateSize[j] = 0;
	}
	for (int i = 0; i < num; i++){
		int temp = indexClis[i];
		if (groupDateSize[temp] < Groupcount){
			copyDouble(points[temp][groupDateSize[temp]], datas[i], colum2);
			groupDateSize[temp]++;
		}
		bool flag = true;
		for (int j = 0; j < size; j++){
			if (groupDateSize[j] < Groupcount){
				flag = false;
				break;
			}
		}
		if (flag){
			return;
		}
	}
}
void getMeanTemp(int *ind, double **datas, double **arr, int num, int colum, int k, int *size){

	int i, j, index;
	for (i = 0; i < k; i++){
		size[i] = 0;
		for (j = 0; j < colum; j++){
			arr[i][j] = 0;
		}
	}
	for (i = 0; i < num; i++){
		index = ind[i];
		for (j = 0; j < colum; j++){
			arr[index][j] += datas[i][j];
		}
		size[index]++;
	}
	for (i = 0; i < k; i++){
		if (size[i]){
			for (j = 0; j < colum; j++){
				arr[i][j] /= size[i];
			}
		}
	}
	return;
}
void getMeanGroup(int *ind, double **datas, double ***arr, double ***point, int num, int colum, int k, int **size, int groupCount, int*indTemp){
	int i, j, index;
	for (i = 0; i < groupCount; i++){
		setZeroInteger(size[i], k);
		for (j = 0; j < k; j++){
			setZeroDouble(arr[i][j], colum);
		}
	}
	//for (i = 0; i < num; i++){
	//	index = ind[i];
	//	for (j = 0; j < colum; j++){
	//		arr[index][j] += datas[i][j];
	//	}
	//	size[index]++;
	//}
	for (i = 0; i < num; i++){
		index = ind[i];
		arrAdd(arr[indTemp[i]][index%k], datas[i], colum);
		size[indTemp[i]][index%k]++;
	}

	for (i = 0; i < groupCount; i++){
		for (j = 0; j < k; j++){
			if (size[i][j]){
				arrDivide(arr[i][j], size[i][j], colum);
			}
			else{
				copyDouble(arr[i][j], point[i][j], colum);
			}
		}
	}

	return;
}
double distEclud(double **data, double **C, int data_i, int C_j, int dataDim) {
	double dist = 0;

	for (int l = 0; l < dataDim; l++) {
		dist += (data[data_i][l] - C[C_j][l]) * (data[data_i][l] - C[C_j][l]);
	}

	return dist;
}
void sccodeInit(double **data, double **C_init, int dataNum, int cluNum, int dataDim) {
	cout << "源代码初始化方法" << endl;
	int minFstDimIdx = rand() % dataNum;
	cout << minFstDimIdx << ":";
	for (int j = 0; j < dataDim; j++) {
		C_init[0][j] = data[minFstDimIdx][j];
	}
	vector<int> nodes(dataNum);
	for (int i = 0; i < dataNum; i++) {
		nodes[i] = i;
	}
	random_shuffle(nodes.begin(), nodes.end());
	//for (int i = 0; i < dataNum;i++){cout<<nodes[i]<<" ";}
	//cout<<endl;
	vector<bool> isin(dataNum, false);
	isin[minFstDimIdx] = true;
	int L = cluNum * 30;
	int nodescount = min(L, dataNum);
	vector<double>sumDisFrom_K_1_PreviousCenter(nodescount, 0.0);

	for (int i = 1; i < cluNum; i++) {
		int nextnode = -1;
		double nextweight = DBL_MIN;

		for (int z = 0; z < nodescount; z++) {
			int node = nodes[z];

			if (isin[node]) {
				continue;
			}  //{ L++;  }

			double newdisFromLastCenter = distEclud(data, C_init, node, i - 1, dataDim);
			//double newdisFromLastCenter = inst->Dis(node, k_1_selectedCenters[i - 1]);
			sumDisFrom_K_1_PreviousCenter[z] = sumDisFrom_K_1_PreviousCenter[z] + newdisFromLastCenter;
			double diff = abs(newdisFromLastCenter - sumDisFrom_K_1_PreviousCenter[z] / (double)i);

			if (diff == 0) {
				diff = 1;
			}

			double curweight = sumDisFrom_K_1_PreviousCenter[z] / diff;

			if (curweight > nextweight) {
				nextweight = curweight;
				nextnode = node;
			}
		}
		for (int j = 0; j < dataDim; j++) {
			C_init[i][j] = data[nextnode][j];
		}
		isin[nextnode] = true;
	}
}


void datasgetKmeans(int colum1, int colum2, int num, int k, int &count, int flag, double**datas1, double***pp, double &count1, int tan, double &localCount, double &noLocalCount){
	try{
		double **point = getArr(k, colum2);
		sccodeInit(datas1, point, num, k, colum2);
		int dataGroupCount = (num - 1) / 5 + 1;
		//逐次迭代
		//存放地i个数据属于第几个中心
		int * listclusterIndex = new int[num];
		int * listclusterIndexNew = new int[k];
		int * listclusterIndexGroupData = new int[num];
		double *means, *p;
		int i, j, m;
		//存放中心更新前每个中心数据个数
		int *size1 = new int[k];
		//存放中心更新后每个中心数据个数
		int *size2 = new int[k];
		//对存放中心个数初始化
		for (i = 0; i < k; i++){
			size1[i] = 0;
		}
		//存放每次更新每类数据变化值
		double **change = new double*[k];
		*change = new double[k*colum2];
		for (i = 0; i < k; i++)
		{
			change[i] = *change + i*colum2;
		}
		double dis, d, secondDis;
		bool f1;
		int *size = new int[k];
		double **iniPoints = new double*[k];
		*iniPoints = new double[k*colum2];
		for (i = 0; i < k; i++)
		{
			iniPoints[i] = *iniPoints + i*colum2;
		}
		double **arr = new double*[k];
		*arr = new double[k*colum2];
		for (i = 0; i < k; i++)
		{
			arr[i] = *arr + i*colum2;
		}
		if (0 == flag){
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = point[i][j];
				}
			}
		}
		else{
			//TODO：随机中心
		}
		int count2 = 0;
		int count3 = 0;
		int *rands = new int[dataGroupCount];
		//int countTemp = 0;
		//double countTemp1 = 0;
		//double countLocal = 0;
		//double countNoLocal = 0;
		//getRamdArr(datas, datas1, rands, dataGroupCount, colum2, num);
		//double *dblArr = new double[dataGroupCount];
		double high = 6;
		int du = pow(dataGroupCount, 1.0 / (high - 1));
		dataGroupCount = pow(du, high - 1);
		int *datasIndex = new int[num];
		double **datas = getArr(dataGroupCount, colum2);
		double **dataArr = getArr(dataGroupCount, colum2);
		int *groupDataSize = new int[dataGroupCount];
		int *listclusterIndexGroupDataTemp = new int[num];
		double ***datasPoints = new double**[dataGroupCount];
		for (i = 0; i < dataGroupCount; i++){
			datasPoints[i] = getArr(du, colum2);
		}
		double ***datasArr = new double**[dataGroupCount];
		for (i = 0; i < dataGroupCount; i++){
			datasArr[i] = getArr(du, colum2);
		}
		setZeroInteger(listclusterIndexGroupDataTemp, num);
		int**groupArrSize = getIntegerArr(dataGroupCount, du);
		int tempIndex = du;
		int tempIndex2 = 1;
		for (int ii = 0; ii < high - 1; ii++){
			count3 = 0;
			getRamdArr(datasPoints, datas1, du, colum2, num, groupDataSize, listclusterIndexGroupDataTemp, tempIndex2);
			while (count3 < 5){
				for (i = 0; i < num; i++){
					dis = DBL_MAX;
					for (j = 0; j < du; j++){
						d = getDistance1(datas1[i], datasPoints[listclusterIndexGroupDataTemp[i]][j], colum2);
						count1++;
						if (d < dis){
							dis = d;
							listclusterIndexGroupData[i] = listclusterIndexGroupDataTemp[i] * du + j;
						}
					}
				}
				if (ii != high - 2){
					getMeanGroup(listclusterIndexGroupData, datas1, datasArr, datasPoints, num, colum2, du, groupArrSize, tempIndex2, listclusterIndexGroupDataTemp);
					for (i = 0; i < tempIndex2; i++){
						for (j = 0; j < du; j++){
							copyDouble(datasPoints[i][j], datasArr[i][j], colum2);
						}
					}
				}
				/*	getMeanTwo(listclusterIndexGroupData, datas1, dataArr, datas, num, colum2, dataGroupCount, dataSize1);
				for (i = 0; i < dataGroupCount; i++){
				for (j = 0; j < colum2; j++){
				datas[i][j] = dataArr[i][j];
				}
				}*/
				count3++;
			}
			copyInteger(listclusterIndexGroupDataTemp, listclusterIndexGroupData, num);
			tempIndex = tempIndex*du;
			tempIndex2 = tempIndex2*du;
		}
		getMeanTemp(listclusterIndexGroupDataTemp, datas1, datas, num, colum2, dataGroupCount, groupDataSize);
		int sums = 0;
		for (int i = 0; i<dataGroupCount; i++){
			sums += groupDataSize[i];

		}
		cout << sums << endl;
		cout << dataGroupCount;
		num = dataGroupCount;
		//中心分组
		int Groupcount = (k - 1) / tan + 1;
		double **groupPoints = new double*[Groupcount];
		*groupPoints = new double[Groupcount*colum2];
		for (i = 0; i < Groupcount; i++)
		{
			groupPoints[i] = *groupPoints + i*colum2;
		}
		srand(time(NULL));
		int RandNum;
		int *arrRand = new int[Groupcount];
		int front = 0;
		int  flag1 = 0, t1 = 0;
		while (1)
		{
			flag1 = 0;
			if (t1 == Groupcount)
				break;

			RandNum = (rand() % k);
			for (i = 0; i < front; i++)
			{
				if (arrRand[i] == RandNum)
					flag1 = 1;
			}
			if (flag1 != 1)
			{
				arrRand[front++] = RandNum;
				for (j = 0; j < colum2; j++){
					groupPoints[t1][j] = iniPoints[RandNum][j];
				}
				t1++;
			}
		}
		int *groupSize = new int[Groupcount];
		for (i = 0; i < Groupcount; i++){
			groupSize[i] = 0;
		}
		int **groupPointsData = new int*[Groupcount];
		*groupPointsData = new int[Groupcount*k];
		for (i = 0; i < Groupcount; i++)
		{
			groupPointsData[i] = *groupPointsData + i*k;
		}
		while (count3 < 4){
			for (i = 0; i < Groupcount; i++){
				size[i] = 0;
			}
			for (i = 0; i < k; i++){
				dis = DBL_MAX;
				for (j = 0; j < Groupcount; j++){
					d = getDistance1(iniPoints[i], groupPoints[j], colum2);
					count1++;
					if (d < dis){
						dis = d;
						listclusterIndexNew[i] = j;
					}
				}
			}
			getMeanTwo(listclusterIndexNew, iniPoints, arr, groupPoints, k, colum2, Groupcount, size);
			for (i = 0; i < Groupcount; i++){
				for (j = 0; j < colum2; j++){
					groupPoints[i][j] = arr[i][j];
				}
			}
			count3++;
		}
		for (i = 0; i < Groupcount; i++){
			size[i] = 0;
		}
		for (i = 0; i < k; i++){
			dis = DBL_MAX;
			for (j = 0; j < Groupcount; j++){
				d = getDistance1(iniPoints[i], groupPoints[j], colum2);
				count1++;
				if (d < dis){
					dis = d;
					listclusterIndexNew[i] = j;
				}
			}
			groupPointsData[listclusterIndexNew[i]][groupSize[listclusterIndexNew[i]]++] = i;
		}
		int m1 = 0;
		for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
				for (int i1 = 0; i1 < colum2; i1++){
					iniPoints[m1][i1] = point[groupPointsData[i][j]][i1];
				}
				m1++;
			}
		}
		m1 = 0;
		for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
				groupPointsData[i][j] = m1;
				m1++;
			}
		}
		bubble_sort(listclusterIndexNew, k);
		while (false){
			for (i = 0; i < k; i++){
				size[i] = 0;
			}
			for (i = 0; i < num; i++){
				dis = DBL_MAX;
				for (j = 0; j < k; j++){
					d = getDistance1(iniPoints[j], datas[i], colum2);
					count1++;
					if (d < dis){
						dis = d;
						listclusterIndex[i] = j;
					}
				}
			}
			f1 = true;
			getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
			for (i = 0; i < k; i++){
				if (!equalsArr(arr[i], iniPoints[i], colum2)){
					f1 = false;
					break;
				}
			}
			count++;
			//与上次相同停止循环
			if (f1 || (count == max)){
				*pp = iniPoints;
				delete[] * arr;
				delete[] arr;
				//delete[] * iniPoints;
				delete[] iniPoints;
				delete[] size;
				return;
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = arr[i][j];
				}
			}
		}
		//最后一组的中心个数
		//int finalGroupCluster = Groupcount == 1 ? k : Groupcount % tan == 0 ? tan : k % tan;
		//存放每个数据对应每组的下届
		double **listsecondminumgroup = new double*[num];
		*listsecondminumgroup = new double[Groupcount*num];
		for (i = 0; i < num; i++)
		{
			listsecondminumgroup[i] = *listsecondminumgroup + i*Groupcount;
		}
		//double *listsecondminum_distance = new double[num];

		for (i = 0; i < k; i++){
			size[i] = 0;
		}

		//存放全局上届 
		double *listsmallest_distance = new double[num];
		for (i = 0; i < num; i++){
			t1 = groupSize[0];
			dis = DBL_MAX;
			secondDis = DBL_MAX;
			listclusterIndex[i] = 0;
			int top = 0, top1 = 1;
			for (j = 0; j < k; j++){
				d = getDistance(iniPoints[j], datas[i], colum2);
				count1++;
				if (d < dis){
					//如果当前点所在的临时中心与第j个中心为同一个群组说明上一次计算的最小值此时为第二小否则说明上一个群组已经跑完该点不会属于上一个群组将上一个次计算的最小值作为上一个群组的下届 
					if (listclusterIndexNew[listclusterIndex[i]] == listclusterIndexNew[j]){
						secondDis = dis;
					}
					else{
						listsecondminumgroup[i][listclusterIndexNew[listclusterIndex[i]]] = dis;
					}
					dis = d;
					listclusterIndex[i] = j;

				}
				//如果d大于dis且小于secondDis说明此时的d为第二小值 
				else if (d < secondDis){
					secondDis = d;
				}
				//如果当前点已经是该群组的最后一个则把第二小的点作为该组的下届如果最终这个点不在这个群组应该把最小值作为下届此时该数据为错误的这会在下一次循环上面代码中重新进行更新 
				if ((j + 1) == t1){
					listsecondminumgroup[i][top++] = secondDis;
					if (top1 != Groupcount){
						t1 = t1 + groupSize[top1++];
					}
					secondDis = DBL_MAX;
				}
			}
			//listsecondminum_distance[i] = secondDis;
			//该数据的上届更新为最小值 
			listsmallest_distance[i] = dis;
			size1[listclusterIndex[i]] += groupDataSize[i];
		}
		//组最大变化距离 
		double groupbiggest_shift = -1;
		//每组的最大变化距离 
		double *listgroupbiggest_shift = new double[Groupcount];
		//每个中心点的变化距离 
		double *listshift = new double[k];

		int top = 0, top1 = 1;
		f1 = true;
		//用传统方法更新中心 
		getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
		int ss = 0;
		for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
				double shift = getDistance(iniPoints[ss], arr[ss], colum2);
				count1++;
				listshift[ss] = shift;
				//找到当群组的最大变化距离 
				groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
				ss++;
			}
			listgroupbiggest_shift[i] = groupbiggest_shift;
			groupbiggest_shift = -1;
		}
		for (i = 0; i < k; i++){
			if (!equalsArr(arr[i], iniPoints[i], colum2)){
				f1 = false;
				break;
			}
		}
		count++;
		//与上次相同停止循环
		if (f1 || (count == max)){
			*pp = arr;
			delete[] listclusterIndex;
			delete[] listclusterIndexNew;
			delete[] size1;
			delete[] size2;
			delete[] * change;
			delete[] change;
			delete[] * iniPoints;
			delete[] iniPoints;
			delete[] * groupPoints;
			delete[] groupPoints;
			delete[] arrRand;
			delete[] groupSize;
			delete[] * listsecondminumgroup;
			delete[] listsecondminumgroup;
			delete[] listsmallest_distance;
			delete[] listgroupbiggest_shift;
			delete[] listshift;
			return;
		}
		for (i = 0; i < k; i++){
			for (j = 0; j < colum2; j++){
				iniPoints[i][j] = arr[i][j];
			}
		}
		/*************************************************/
		//第三步
		/*************************************************/
		double *listsecondminumgroup_old = new double[Groupcount];
		int *listg = new int[Groupcount];
		while (true){
			/*for (i = 0; i < k; i++){
			size[i] = 0;
			}*/
			for (i = 0; i < k; i++){
				size2[i] = size1[i];
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					change[i][j] = 0;
				}
			}
			int g = 0;
			for (i = 0; i < num; i++){
				int top = 0;
				//用以前的组上届加上改点所在中心的变化距离之和作为新的组上届 
				listsmallest_distance[i] = listsmallest_distance[i] + listshift[listclusterIndex[i]];
				for (j = 0; j < Groupcount; j++){
					listsecondminumgroup_old[j] = listsecondminumgroup[i][j];
				}
				for (j = 0; j < Groupcount; j++){
					//用上一次的组下届与最大变化距离之差作为新的组下届 
					listsecondminumgroup[i][j] -= listgroupbiggest_shift[j];
				}
				int listclusterIndex_old = listclusterIndex[i];
				//int indd = listsmallest_distance[i];
				//如果全局组下届已经大于组上届说明改点不需要改变中心跳出循环否则的话需要继续过滤 
				if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
					//将组上届缩小为改点到当前中心的距离
					double ubTight = getDistance(datas[i], iniPoints[listclusterIndex[i]], colum2);
					count1++;
					listsmallest_distance[i] = ubTight;
					//如果全局下届依旧小于缩小后的上届那么要选出每一个小于上届的组中心点有可能存在其中 
					if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
						for (j = 0; j < Groupcount; j++){

							if (listsecondminumgroup[i][j] < (listsmallest_distance[i])){
								int qqq = groupSize[j];
								noLocalCount += (double)qqq;
								//将组下届更新为正无穷 
								listsecondminumgroup[i][j] = DBL_MAX;
								listg[top++] = j;
							}
						}
						flag = listclusterIndex[i];
						for (j = 0; j < top; j++){
							int t = listg[j];
							int A = groupSize[t];

							for (m = 0; m < A; m++){
								//由于当前最小值为该点当前的旧中心距离且在下面会与当前最小值比较如果不大于它不会进行更新所以不必对它进行比较 
								if (groupPointsData[t][m] != listclusterIndex_old){
									/*	dis = DBL_MAX;
									d = getDistance(iniPoints[tan * t + m], datas[i], colum2);
									if (d < dis){
									dis = d;
									listclusterIndex[i] = tan * t + m;
									}*/
									//根据引理二 
									if (listsecondminumgroup[i][t] >listsecondminumgroup_old[t] - listshift[groupPointsData[t][m]]){

										double lbtight = getDistance(datas[i], iniPoints[groupPointsData[t][m]], colum2);
										localCount++;
										if (lbtight < listsmallest_distance[i]){
											//如果组没有变那么该组第二小要变成上一次的最小值因为当前认为数据属于该组下届需要为第二小 
											if (listclusterIndexNew[flag] == t){
												//if (listclusterIndex[i] / tan == (tan * t + m) / tan){
												listsecondminumgroup[i][t] = listsmallest_distance[i];
											}
											else{
												//listsecondminumgroup[i][listclusterIndex[i] / tan] = listsmallest_distance[i];
												//listsecondminumgroup[i][flag / tan] = listsmallest_distance[i];
												listsecondminumgroup[i][listclusterIndexNew[flag]] = listsmallest_distance[i];
											}
											//更新上届 
											listsmallest_distance[i] = lbtight;
											//更新待定中心 
											flag = groupPointsData[t][m];
											//listclusterIndex[i] = (tan * t + m);
										}
										else if (lbtight < listsecondminumgroup[i][t]){
											listsecondminumgroup[i][t] = lbtight;
										}
									}
								}
							}
						}
						//此处要对变化的数据差值进行统计以优化更新中心的算法 
						int &q = listclusterIndex[i];
						if (flag != q){
							size1[flag] += groupDataSize[i];
							size1[q] -= groupDataSize[i];
							for (j = 0; j < colum2; j++){
								change[flag][j] += datas[i][j] * groupDataSize[i];
								change[q][j] -= datas[i][j] * groupDataSize[i];
							}
						}
						q = flag;
					}
				}
			}

			for (i = 0; i < k; i++){
				double *q = arr[i];
				double *p = change[i];
				int s1 = size1[i];
				int s2 = size2[i];
				for (j = 0; j < colum2; j++){
					double &t = q[j];
					if (s1 != 0)
					{
						t = (t*s2 + p[j]) / s1;
					}

				}
			}
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/*	if (count == 3){
			double **iniTemp = getArr(k, colum2);
			double ** arrTemp = getArr(k, colum2);

			int front = 0;
			int  flag1 = 0, t1 = 0;
			while (1)
			{
			flag1 = 0;
			if (t1 == Groupcount)
			break;

			RandNum = (rand() % k);
			for (i = 0; i < front; i++)
			{
			if (arrRand[i] == RandNum)
			flag1 = 1;
			}
			if (flag1 != 1)
			{
			arrRand[front++] = RandNum;
			for (j = 0; j < colum2; j++){
			groupPoints[t1][j] = arr[RandNum][j];
			}
			t1++;
			}
			}
			int* listindexTemp = new int[num];
			int *groupSizeTemp = new int[Groupcount];
			for (i = 0; i < Groupcount; i++){
			groupSizeTemp[i] = 0;
			}
			int **groupPointsDataTemp = new int*[Groupcount];
			*groupPointsDataTemp = new int[Groupcount*k];
			for (i = 0; i < Groupcount; i++)
			{
			groupPointsDataTemp[i] = *groupPointsDataTemp + i*k;
			}
			int * listclusterIndexNewTemp = new int[k];
			double** secondeTemp = getArr(num, Groupcount);
			setZeroInteger(groupSize, Groupcount);
			count3 = 0;
			while (count3 < 4){
			for (i = 0; i < Groupcount; i++){
			size[i] = 0;
			}
			for (i = 0; i < k; i++){
			dis = DBL_MAX;
			for (j = 0; j < Groupcount; j++){
			d = getDistance1(arr[i], groupPoints[j], colum2);
			count1++;
			if (d < dis){
			dis = d;
			listclusterIndexNewTemp[i] = j;
			}
			}
			}
			getMeanTwo(listclusterIndexNewTemp, arr, arrTemp, groupPoints, k, colum2, Groupcount, size);
			for (i = 0; i < Groupcount; i++){
			for (j = 0; j < colum2; j++){
			groupPoints[i][j] = arrTemp[i][j];
			}
			}
			count3++;
			}
			for (i = 0; i < Groupcount; i++){
			size[i] = 0;
			}
			for (i = 0; i < k; i++){
			dis = DBL_MAX;
			for (j = 0; j < Groupcount; j++){
			d = getDistance1(arr[i], groupPoints[j], colum2);
			count1++;
			if (d < dis){
			dis = d;
			listclusterIndexNewTemp[i] = j;
			}
			}
			groupPointsData[listclusterIndexNewTemp[i]][groupSize[listclusterIndexNewTemp[i]]++] = i;
			}

			for (int ii = 0; ii < num; ii++){
			for (int jj = 0; jj < Groupcount; jj++){
			secondeTemp[ii][jj] = DBL_MAX;
			}
			}
			for (int ii = 0; ii < num; ii++){
			for (int jj = 0; jj < k; jj++){
			if (secondeTemp[ii][listclusterIndexNewTemp[jj]] > listsecondminumgroup[ii][listclusterIndexNew[jj]]){
			secondeTemp[ii][listclusterIndexNewTemp[jj]] = listsecondminumgroup[ii][listclusterIndexNew[jj]];
			}
			}
			}
			for (int i = 0; i < k; i++){
			copyDouble(iniTemp[i], iniPoints[i], colum2);
			copyDouble(arrTemp[i], arr[i], colum2);
			}
			int *sizeTemp = new int[k];
			int *arrIndex = new int[num];
			int m1 = 0;
			for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
			for (int i1 = 0; i1 < colum2; i1++){
			iniPoints[m1][i1] = iniTemp[groupPointsData[i][j]][i1];
			arr[m1][i1] = arrTemp[groupPointsData[i][j]][i1];
			}
			sizeTemp[m1] = size1[groupPointsData[i][j]];
			arrIndex[groupPointsData[i][j]] = m1;
			m1++;
			}
			}
			for (i = 0; i < num; i++){
			listindexTemp[i] = arrIndex[listclusterIndex[i]];
			}
			m1 = 0;
			for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
			groupPointsData[i][j] = m1;
			m1++;
			}
			}
			bubble_sort(listclusterIndexNewTemp, k);
			copyInteger(listclusterIndex, listindexTemp, num);
			copyInteger(listclusterIndexNew, listclusterIndexNewTemp, k);
			copyInteger(size1, sizeTemp, k);
			for (i = 0; i < num; i++){
			copyDouble(secondeTemp[i], listsecondminumgroup[i], Groupcount);
			}

			}*/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			/********************************************************************************************************************************/
			top = 0, top1 = 1;
			f1 = true;
			int ss = 0;
			for (i = 0; i < Groupcount; i++){
				for (j = 0; j < groupSize[i]; j++){
					double shift = getDistance(iniPoints[ss], arr[ss], colum2);
					count1++;
					listshift[ss] = shift;
					//找到当群组的最大变化距离 
					groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
					ss++;
				}
				listgroupbiggest_shift[i] = groupbiggest_shift;
				groupbiggest_shift = -1;
			}

			for (i = 0; i < k; i++){
				if (!equalsArr(arr[i], iniPoints[i], colum2)){
					f1 = false;
					break;
				}
			}
			count++;
			//与上次相同停止循环
			if (f1 || (count == max)){
				*pp = arr;
				delete[] listclusterIndex;
				delete[] listclusterIndexNew;
				delete[] size1;
				delete[] size2;
				delete[] * change;
				delete[] change;
				delete[] * iniPoints;
				delete[] iniPoints;
				delete[] * groupPoints;
				delete[] groupPoints;
				delete[] arrRand;
				delete[] groupSize;
				delete[] * listsecondminumgroup;
				delete[] listsecondminumgroup;
				delete[] listsmallest_distance;
				delete[] listgroupbiggest_shift;
				delete[] listshift;
				delete[] listsecondminumgroup_old;
				delete[] listg;
				return;
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = arr[i][j];
				}
			}
		}
	}
	catch (char *str)
	{
		cout << str << endl;
	}
	catch (int i)
	{
		cout << i << endl;
	}
	/*************************************************/
}
double calObjFuncForGivenCents(double **data, double **V, int data_num, int data_dim, int clu_num)
{
	double objFunc = 0;
	double tmpMinDistSquare; // 距离平方
	double tmpDistSquare;  // 距离平方
	for (int i = 0; i < data_num; i++)
	{
		tmpMinDistSquare = DBL_MAX;
		for (int j = 0; j < clu_num; j++)
		{
			tmpDistSquare = 0;
			for (int l = 0; l < data_dim; l++)
			{
				tmpDistSquare += ((data[i][l] - V[j][l])*(data[i][l] - V[j][l]));
			}
			if (tmpDistSquare < tmpMinDistSquare)
			{
				tmpMinDistSquare = tmpDistSquare;
			}
		}
		objFunc += tmpMinDistSquare;
	}

	return objFunc;
}
double min(double *arr, int length){
	double t = arr[0];
	for (int i = 1; i < length; i++){
		if (t > arr[i]){
			t = arr[i];
		}
	}
	return t;
}
void getMeanTwo(int *ind, double **datas, double **arr, double**point, int num, int colum, int k, int *size){

	int i, j, index;
	for (i = 0; i < k; i++){
		size[i] = 0;
		for (j = 0; j < colum; j++){
			arr[i][j] = 0;
		}
	}
	for (i = 0; i < num; i++){
		index = ind[i];
		for (j = 0; j < colum; j++){
			arr[index][j] += datas[i][j];
		}
		size[index]++;
	}
	for (i = 0; i < k; i++){
		if (size[i]){
			for (j = 0; j < colum; j++){
				arr[i][j] /= size[i];
			}
		}
		else{
			for (j = 0; j < colum; j++){
				arr[i][j] = point[i][j];
			}
		}

	}
	return;
}
inline bool equalsArr(double *vec1, double *vec2, int colum){
	int i;
	for (i = 0; i < colum; i++){
		if (vec1[i] != vec2[i]){
			return false;
		}
	}
	return true;
}
void getArr(int a, int b, double **c){
	*c = new double[a*b];
	for (int i = 0; i < a; i++)
	{
		c[i] = *c + i*b;
	}
}

inline void getMean(double**datas, int *data1, double data2[], int size, int colum){
	int i, j;
	for (i = 0; i < colum; i++){
		double sum = 0.0;
		for (j = 0; j < size; j++){
			sum += datas[data1[j]][i];
		}
		data2[i] = sum / size;
	}
}
void datasgetKmeansYinYang(int colum1, int colum2, int num, int k, int &count, double **point, int flag, double**datas, double***pp, double &count1, int tan, double &localCount, double &noLocalCount, int *listclusterIndex){
	try{

		//逐次迭代
		//存放地i个数据属于第几个中心
		int * listclusterIndexNew = new int[k];
		double *means, *p;
		int i, j, m;
		//存放中心更新前每个中心数据个数
		int *size1 = new int[k];
		//存放中心更新后每个中心数据个数
		int *size2 = new int[k];
		//对存放中心个数初始化
		for (i = 0; i < k; i++){
			size1[i] = 0;
		}
		//存放每次更新每类数据变化值
		double **change = new double*[k];
		*change = new double[k*colum2];
		for (i = 0; i<k; i++)
		{
			change[i] = *change + i*colum2;
		}
		double dis, d, secondDis;
		bool f1;
		int *size = new int[k];
		double **iniPoints = new double*[k];
		*iniPoints = new double[k*colum2];
		for (i = 0; i<k; i++)
		{
			iniPoints[i] = *iniPoints + i*colum2;
		}

		double **arr = new double*[k];
		*arr = new double[k*colum2];
		for (i = 0; i<k; i++)
		{
			arr[i] = *arr + i*colum2;
		}
		if (0 == flag){
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = point[i][j];
				}
			}
		}
		else{
			//TODO：随机中心
		}
		int count2 = 0;
		int count3 = 0;





		//中心分组
		int Groupcount = (k - 1) / tan + 1;

		double **groupPoints = new double*[Groupcount];
		*groupPoints = new double[Groupcount*colum2];
		for (i = 0; i<Groupcount; i++)
		{
			groupPoints[i] = *groupPoints + i*colum2;
		}
		srand(time(NULL));
		int RandNum;
		int *arrRand = new int[Groupcount];
		int front = 0;
		int  flag1 = 0, t1 = 0;
		while (1)
		{
			flag1 = 0;
			if (t1 == Groupcount)
				break;

			RandNum = (rand() % k);
			for (i = 0; i < front; i++)
			{
				if (arrRand[i] == RandNum)
					flag1 = 1;
			}
			if (flag1 != 1)
			{
				arrRand[front++] = RandNum;
				for (j = 0; j<colum2; j++){
					groupPoints[t1][j] = iniPoints[RandNum][j];
				}
				t1++;
			}
		}
		int *groupSize = new int[Groupcount];
		for (i = 0; i < Groupcount; i++){
			groupSize[i] = 0;
		}
		int **groupPointsData = new int*[Groupcount];
		*groupPointsData = new int[Groupcount*k];
		for (i = 0; i<Groupcount; i++)
		{
			groupPointsData[i] = *groupPointsData + i*k;
		}
		while (count3<4){
			for (i = 0; i < Groupcount; i++){
				size[i] = 0;
			}
			for (i = 0; i < k; i++){
				dis = DBL_MAX;
				for (j = 0; j < Groupcount; j++){
					d = getDistance1(iniPoints[i], groupPoints[j], colum2);
					count1++;
					if (d < dis){
						dis = d;
						listclusterIndexNew[i] = j;
					}
				}
			}
			getMeanTwo(listclusterIndexNew, iniPoints, arr, groupPoints, k, colum2, Groupcount, size);
			for (i = 0; i < Groupcount; i++){
				for (j = 0; j<colum2; j++){
					groupPoints[i][j] = arr[i][j];
				}
			}
			count3++;
		}
		for (i = 0; i < Groupcount; i++){
			size[i] = 0;
		}
		for (i = 0; i < k; i++){
			dis = DBL_MAX;
			for (j = 0; j < Groupcount; j++){
				d = getDistance1(iniPoints[i], groupPoints[j], colum2);
				count1++;
				if (d < dis){
					dis = d;
					listclusterIndexNew[i] = j;
				}
			}
			groupPointsData[listclusterIndexNew[i]][groupSize[listclusterIndexNew[i]]++] = i;
		}
		int m1 = 0;
		for (i = 0; i<Groupcount; i++){
			for (j = 0; j<groupSize[i]; j++){
				for (int i1 = 0; i1 < colum2; i1++){
					iniPoints[m1][i1] = point[groupPointsData[i][j]][i1];
				}
				m1++;
			}
		}
		m1 = 0;
		for (i = 0; i<Groupcount; i++){
			for (j = 0; j<groupSize[i]; j++){
				groupPointsData[i][j] = m1;
				m1++;
			}
		}
		bubble_sort(listclusterIndexNew, k);
		while (false){
			for (i = 0; i < k; i++){
				size[i] = 0;
			}
			for (i = 0; i < num; i++){
				dis = DBL_MAX;
				for (j = 0; j < k; j++){
					d = getDistance1(iniPoints[j], datas[i], colum2);
					count1++;
					if (d < dis){
						dis = d;
						listclusterIndex[i] = j;
					}
				}
			}
			f1 = true;
			getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
			for (i = 0; i < k; i++){
				if (!equalsArr(arr[i], iniPoints[i], colum2)){
					f1 = false;
					break;
				}
			}
			count++;
			//与上次相同停止循环
			if (f1 || (count == max)){
				*pp = iniPoints;
				delete[] * arr;
				delete[] arr;
				//delete[] * iniPoints;
				delete[] iniPoints;
				delete[] size;
				return;
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = arr[i][j];
				}
			}
		}
		//最后一组的中心个数
		//int finalGroupCluster = Groupcount == 1 ? k : Groupcount % tan == 0 ? tan : k % tan;
		//存放每个数据对应每组的下届
		double **listsecondminumgroup = new double*[num];
		*listsecondminumgroup = new double[Groupcount*num];
		for (i = 0; i<num; i++)
		{
			listsecondminumgroup[i] = *listsecondminumgroup + i*Groupcount;
		}
		//double *listsecondminum_distance = new double[num];

		for (i = 0; i < k; i++){
			size[i] = 0;
		}

		//存放全局上届 
		double *listsmallest_distance = new double[num];
		for (i = 0; i < num; i++){
			t1 = groupSize[0];
			dis = DBL_MAX;
			secondDis = DBL_MAX;
			listclusterIndex[i] = 0;
			int top = 0, top1 = 1;
			for (j = 0; j < k; j++){
				d = getDistance(iniPoints[j], datas[i], colum2);
				count1++;
				if (d < dis){
					//如果当前点所在的临时中心与第j个中心为同一个群组说明上一次计算的最小值此时为第二小否则说明上一个群组已经跑完该点不会属于上一个群组将上一个次计算的最小值作为上一个群组的下届 
					if (listclusterIndexNew[listclusterIndex[i]] == listclusterIndexNew[j]){
						secondDis = dis;
					}
					else{
						listsecondminumgroup[i][listclusterIndexNew[listclusterIndex[i]]] = dis;
					}
					dis = d;
					listclusterIndex[i] = j;

				}
				//如果d大于dis且小于secondDis说明此时的d为第二小值 
				else if (d<secondDis){
					secondDis = d;
				}
				//如果当前点已经是该群组的最后一个则把第二小的点作为该组的下届如果最终这个点不在这个群组应该把最小值作为下届此时该数据为错误的这会在下一次循环上面代码中重新进行更新 
				if ((j + 1) == t1){
					listsecondminumgroup[i][top++] = secondDis;
					if (top1 != Groupcount){
						t1 = t1 + groupSize[top1++];
					}
					secondDis = DBL_MAX;
				}
			}
			//listsecondminum_distance[i] = secondDis;
			//该数据的上届更新为最小值 
			listsmallest_distance[i] = dis;
			size1[listclusterIndex[i]]++;
		}
		//组最大变化距离 
		double groupbiggest_shift = -1;
		//每组的最大变化距离 
		double *listgroupbiggest_shift = new double[Groupcount];
		//每个中心点的变化距离 
		double *listshift = new double[k];

		int top = 0, top1 = 1;
		f1 = true;
		//用传统方法更新中心 
		getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
		t1 = groupSize[0];
		for (i = 0; i < k; i++){
			//计算中心点变化距离 
			double shift = getDistance(iniPoints[i], arr[i], colum2);
			count1++;
			listshift[i] = shift;
			//找到当群组的最大变化距离 
			groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
			//如果到该组最后一个中心点则要更新组最大距离并且把最大变化距离初始化 
			if ((i + 1) == t1){
				listgroupbiggest_shift[top++] = groupbiggest_shift;
				groupbiggest_shift = -1;
				if (top1 != Groupcount){
					t1 = t1 + groupSize[top1++];
				}

			}
		}
		for (i = 0; i < k; i++){
			if (!equalsArr(arr[i], iniPoints[i], colum2)){
				f1 = false;
				break;
			}
		}
		count++;
		//与上次相同停止循环
		if (f1 || (count == max)){
			*pp = arr;
			delete[] listclusterIndex;
			delete[] listclusterIndexNew;
			delete[] size1;
			delete[] size2;
			delete[] * change;
			delete[] change;
			delete[] * iniPoints;
			delete[] iniPoints;
			delete[] * groupPoints;
			delete[] groupPoints;
			delete[] arrRand;
			delete[] groupSize;
			delete[] * listsecondminumgroup;
			delete[] listsecondminumgroup;
			delete[] listsmallest_distance;
			delete[] listgroupbiggest_shift;
			delete[] listshift;
			return;
		}
		for (i = 0; i < k; i++){
			for (j = 0; j < colum2; j++){
				iniPoints[i][j] = arr[i][j];
			}
		}
		/*************************************************/
		//第三步
		/*************************************************/
		double *listsecondminumgroup_old = new double[Groupcount];
		int *listg = new int[Groupcount];
		while (true){
			/*for (i = 0; i < k; i++){
			size[i] = 0;
			}*/
			for (i = 0; i < k; i++){
				size2[i] = size1[i];
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					change[i][j] = 0;
				}
			}
			int g = 0;
			for (i = 0; i < num; i++){
				int top = 0;
				//用以前的组上届加上改点所在中心的变化距离之和作为新的组上届 
				listsmallest_distance[i] = listsmallest_distance[i] + listshift[listclusterIndex[i]];
				for (j = 0; j < Groupcount; j++){
					listsecondminumgroup_old[j] = listsecondminumgroup[i][j];
				}
				for (j = 0; j < Groupcount; j++){
					//用上一次的组下届与最大变化距离之差作为新的组下届 
					listsecondminumgroup[i][j] -= listgroupbiggest_shift[j];
				}
				int listclusterIndex_old = listclusterIndex[i];
				//int indd = listsmallest_distance[i];
				//如果全局组下届已经大于组上届说明改点不需要改变中心跳出循环否则的话需要继续过滤 
				if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
					//将组上届缩小为改点到当前中心的距离
					double ubTight = getDistance(datas[i], iniPoints[listclusterIndex[i]], colum2);
					count1++;
					listsmallest_distance[i] = ubTight;
					//如果全局下届依旧小于缩小后的上届那么要选出每一个小于上届的组中心点有可能存在其中 
					if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
						for (j = 0; j < Groupcount; j++){

							if (listsecondminumgroup[i][j] < (listsmallest_distance[i])){
								int qqq = groupSize[j];
								noLocalCount += (double)qqq;
								//将组下届更新为正无穷 
								listsecondminumgroup[i][j] = DBL_MAX;
								listg[top++] = j;
							}
						}
						flag = listclusterIndex[i];
						for (j = 0; j < top; j++){
							int t = listg[j];
							int A = groupSize[t];

							for (m = 0; m < A; m++){
								//由于当前最小值为该点当前的旧中心距离且在下面会与当前最小值比较如果不大于它不会进行更新所以不必对它进行比较 
								if (groupPointsData[t][m] != listclusterIndex_old){
									/*	dis = DBL_MAX;
									d = getDistance(iniPoints[tan * t + m], datas[i], colum2);
									if (d < dis){
									dis = d;
									listclusterIndex[i] = tan * t + m;
									}*/
									//根据引理二 
									if (listsecondminumgroup[i][t] >listsecondminumgroup_old[t] - listshift[groupPointsData[t][m]]){

										double lbtight = getDistance(datas[i], iniPoints[groupPointsData[t][m]], colum2);
										localCount++;
										if (lbtight < listsmallest_distance[i]){
											//如果组没有变那么该组第二小要变成上一次的最小值因为当前认为数据属于该组下届需要为第二小 
											if (listclusterIndexNew[flag] == t){
												//if (listclusterIndex[i] / tan == (tan * t + m) / tan){
												listsecondminumgroup[i][t] = listsmallest_distance[i];
											}
											else{
												//listsecondminumgroup[i][listclusterIndex[i] / tan] = listsmallest_distance[i];
												//listsecondminumgroup[i][flag / tan] = listsmallest_distance[i];
												listsecondminumgroup[i][listclusterIndexNew[flag]] = listsmallest_distance[i];
											}
											//更新上届 
											listsmallest_distance[i] = lbtight;
											//更新待定中心 
											flag = groupPointsData[t][m];
											//listclusterIndex[i] = (tan * t + m);
										}
										else if (lbtight < listsecondminumgroup[i][t]){
											listsecondminumgroup[i][t] = lbtight;
										}
									}
								}
							}
						}
						//此处要对变化的数据差值进行统计以优化更新中心的算法 
						int &q = listclusterIndex[i];
						if (flag != q){
							size1[flag]++;
							size1[q]--;
							for (j = 0; j < colum2; j++){
								change[flag][j] += datas[i][j];
								change[q][j] -= datas[i][j];
							}
						}
						q = flag;
					}
				}
			}

			for (i = 0; i < k; i++){
				double *q = arr[i];
				double *p = change[i];
				int s1 = size1[i];
				int s2 = size2[i];
				for (j = 0; j < colum2; j++){
					double &t = q[j];
					if (s1 != 0)
					{
						t = (t*s2 + p[j]) / s1;
					}

				}
			}
			top = 0, top1 = 1;
			f1 = true;
			int t1 = groupSize[0];
			for (i = 0; i < k; i++){
				//计算中心点变化距离 
				double shift = getDistance(iniPoints[i], arr[i], colum2);
				count1++;
				listshift[i] = shift;
				//找到当群组的最大变化距离 
				groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
				//如果到该组最后一个中心点则要更新组最大距离并且把最大变化距离初始化 
				if ((i + 1) == t1){
					listgroupbiggest_shift[top++] = groupbiggest_shift;
					groupbiggest_shift = -1;
					if (top1 != Groupcount){
						t1 = t1 + groupSize[top1++];
					}
				}
			}

			for (i = 0; i < k; i++){
				if (!equalsArr(arr[i], iniPoints[i], colum2)){
					f1 = false;
					break;
				}
			}
			count++;
			//与上次相同停止循环
			if (f1 || (count == 5)){
				*pp = arr;
				/*	delete[] listclusterIndex;*/
				delete[] listclusterIndexNew;
				delete[] size1;
				delete[] size2;
				delete[] * change;
				delete[] change;
				delete[] * iniPoints;
				delete[] iniPoints;
				delete[] * groupPoints;
				delete[] groupPoints;
				delete[] arrRand;
				delete[] groupSize;
				delete[] * listsecondminumgroup;
				delete[] listsecondminumgroup;
				delete[] listsmallest_distance;
				delete[] listgroupbiggest_shift;
				delete[] listshift;
				delete[] listsecondminumgroup_old;
				delete[] listg;
				return;
			}
			for (i = 0; i < k; i++){
				for (j = 0; j < colum2; j++){
					iniPoints[i][j] = arr[i][j];
				}
			}
		}
	}
	catch (char *str)
	{
		cout << str << endl;
	}
	catch (int i)
	{
		cout << i << endl;
	}
	/*************************************************/
}




