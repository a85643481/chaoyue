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
	if (n < 0)// ������
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
	cout << "�����������ļ�����" << endl;
	cin >> dataFile;
	cout << "���������ݸ���" << endl;
	cin >> num;
	cout << "����������ά��" << endl;
	cin >> dim;
	cout << "������������k" << endl;
	cin >> K;
	cout << "���������д���" << endl;
	cin >> times;
	cout << "������Ⱥ�����k/x�е�xֵ" << endl;
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
	string path = dataFile + "_c=" + to_String(K) + "_dim=" + to_String(dim) + "_dataNum=" + to_String(num) + "_maxRunTimes=" + to_String(times) + "_yinyang" + "(" + "t=" + to_String((K - 1) / x + 1) + "=" + to_String(K) + "��" + to_String(x) + ")_" + ch + ".txt";
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
	fout << setw(35) << "      [���پ������������� (����)]";
	fout << setw(35) << "      [���پ������������� (��local)]";
	fout << setw(35) << "             timeAll(����)";
	fout << setw(35) << "                   timePerIter(����)" << endl;;
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
	cout << "�������";
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
//�������ʼ��0
void setZeroDouble(double*a, int n){
	for (int i = 0; i < n; i++){
		a[i] = 0;
	}
}
//�������ʼ��0
void setZeroInteger(int*a, int n){
	for (int i = 0; i < n; i++){
		a[i] = 0;
	}
}
//����һ������
void copyDouble(double*a, double*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//����һ������
void copyInteger(int*a, int*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//����һ������
void copyBool(bool*a, bool*b, int n){
	for (int i = 0; i < n; i++){
		a[i] = b[i];
	}
}
//��������һ������ֵ����һ������
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
//�����������
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
	cout << "Դ�����ʼ������" << endl;
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
		//��ε���
		//��ŵ�i���������ڵڼ�������
		int * listclusterIndex = new int[num];
		int * listclusterIndexNew = new int[k];
		int * listclusterIndexGroupData = new int[num];
		double *means, *p;
		int i, j, m;
		//������ĸ���ǰÿ���������ݸ���
		int *size1 = new int[k];
		//������ĸ��º�ÿ���������ݸ���
		int *size2 = new int[k];
		//�Դ�����ĸ�����ʼ��
		for (i = 0; i < k; i++){
			size1[i] = 0;
		}
		//���ÿ�θ���ÿ�����ݱ仯ֵ
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
			//TODO���������
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
		//���ķ���
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
			//���ϴ���ֹͬͣѭ��
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
		//���һ������ĸ���
		//int finalGroupCluster = Groupcount == 1 ? k : Groupcount % tan == 0 ? tan : k % tan;
		//���ÿ�����ݶ�Ӧÿ����½�
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

		//���ȫ���Ͻ� 
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
					//�����ǰ�����ڵ���ʱ�������j������Ϊͬһ��Ⱥ��˵����һ�μ������Сֵ��ʱΪ�ڶ�С����˵����һ��Ⱥ���Ѿ�����õ㲻��������һ��Ⱥ�齫��һ���μ������Сֵ��Ϊ��һ��Ⱥ����½� 
					if (listclusterIndexNew[listclusterIndex[i]] == listclusterIndexNew[j]){
						secondDis = dis;
					}
					else{
						listsecondminumgroup[i][listclusterIndexNew[listclusterIndex[i]]] = dis;
					}
					dis = d;
					listclusterIndex[i] = j;

				}
				//���d����dis��С��secondDis˵����ʱ��dΪ�ڶ�Сֵ 
				else if (d < secondDis){
					secondDis = d;
				}
				//�����ǰ���Ѿ��Ǹ�Ⱥ������һ����ѵڶ�С�ĵ���Ϊ������½������������㲻�����Ⱥ��Ӧ�ð���Сֵ��Ϊ�½��ʱ������Ϊ������������һ��ѭ��������������½��и��� 
				if ((j + 1) == t1){
					listsecondminumgroup[i][top++] = secondDis;
					if (top1 != Groupcount){
						t1 = t1 + groupSize[top1++];
					}
					secondDis = DBL_MAX;
				}
			}
			//listsecondminum_distance[i] = secondDis;
			//�����ݵ��Ͻ����Ϊ��Сֵ 
			listsmallest_distance[i] = dis;
			size1[listclusterIndex[i]] += groupDataSize[i];
		}
		//�����仯���� 
		double groupbiggest_shift = -1;
		//ÿ������仯���� 
		double *listgroupbiggest_shift = new double[Groupcount];
		//ÿ�����ĵ�ı仯���� 
		double *listshift = new double[k];

		int top = 0, top1 = 1;
		f1 = true;
		//�ô�ͳ������������ 
		getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
		int ss = 0;
		for (i = 0; i < Groupcount; i++){
			for (j = 0; j < groupSize[i]; j++){
				double shift = getDistance(iniPoints[ss], arr[ss], colum2);
				count1++;
				listshift[ss] = shift;
				//�ҵ���Ⱥ������仯���� 
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
		//���ϴ���ֹͬͣѭ��
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
		//������
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
				//����ǰ�����Ͻ���ϸĵ��������ĵı仯����֮����Ϊ�µ����Ͻ� 
				listsmallest_distance[i] = listsmallest_distance[i] + listshift[listclusterIndex[i]];
				for (j = 0; j < Groupcount; j++){
					listsecondminumgroup_old[j] = listsecondminumgroup[i][j];
				}
				for (j = 0; j < Groupcount; j++){
					//����һ�ε����½������仯����֮����Ϊ�µ����½� 
					listsecondminumgroup[i][j] -= listgroupbiggest_shift[j];
				}
				int listclusterIndex_old = listclusterIndex[i];
				//int indd = listsmallest_distance[i];
				//���ȫ�����½��Ѿ��������Ͻ�˵���ĵ㲻��Ҫ�ı���������ѭ������Ļ���Ҫ�������� 
				if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
					//�����Ͻ���СΪ�ĵ㵽��ǰ���ĵľ���
					double ubTight = getDistance(datas[i], iniPoints[listclusterIndex[i]], colum2);
					count1++;
					listsmallest_distance[i] = ubTight;
					//���ȫ���½�����С����С����Ͻ���ôҪѡ��ÿһ��С���Ͻ�������ĵ��п��ܴ������� 
					if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
						for (j = 0; j < Groupcount; j++){

							if (listsecondminumgroup[i][j] < (listsmallest_distance[i])){
								int qqq = groupSize[j];
								noLocalCount += (double)qqq;
								//�����½����Ϊ������ 
								listsecondminumgroup[i][j] = DBL_MAX;
								listg[top++] = j;
							}
						}
						flag = listclusterIndex[i];
						for (j = 0; j < top; j++){
							int t = listg[j];
							int A = groupSize[t];

							for (m = 0; m < A; m++){
								//���ڵ�ǰ��СֵΪ�õ㵱ǰ�ľ����ľ�������������뵱ǰ��Сֵ�Ƚ������������������и������Բ��ض������бȽ� 
								if (groupPointsData[t][m] != listclusterIndex_old){
									/*	dis = DBL_MAX;
									d = getDistance(iniPoints[tan * t + m], datas[i], colum2);
									if (d < dis){
									dis = d;
									listclusterIndex[i] = tan * t + m;
									}*/
									//��������� 
									if (listsecondminumgroup[i][t] >listsecondminumgroup_old[t] - listshift[groupPointsData[t][m]]){

										double lbtight = getDistance(datas[i], iniPoints[groupPointsData[t][m]], colum2);
										localCount++;
										if (lbtight < listsmallest_distance[i]){
											//�����û�б���ô����ڶ�СҪ�����һ�ε���Сֵ��Ϊ��ǰ��Ϊ�������ڸ����½���ҪΪ�ڶ�С 
											if (listclusterIndexNew[flag] == t){
												//if (listclusterIndex[i] / tan == (tan * t + m) / tan){
												listsecondminumgroup[i][t] = listsmallest_distance[i];
											}
											else{
												//listsecondminumgroup[i][listclusterIndex[i] / tan] = listsmallest_distance[i];
												//listsecondminumgroup[i][flag / tan] = listsmallest_distance[i];
												listsecondminumgroup[i][listclusterIndexNew[flag]] = listsmallest_distance[i];
											}
											//�����Ͻ� 
											listsmallest_distance[i] = lbtight;
											//���´������� 
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
						//�˴�Ҫ�Ա仯�����ݲ�ֵ����ͳ�����Ż��������ĵ��㷨 
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
					//�ҵ���Ⱥ������仯���� 
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
			//���ϴ���ֹͬͣѭ��
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
	double tmpMinDistSquare; // ����ƽ��
	double tmpDistSquare;  // ����ƽ��
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

		//��ε���
		//��ŵ�i���������ڵڼ�������
		int * listclusterIndexNew = new int[k];
		double *means, *p;
		int i, j, m;
		//������ĸ���ǰÿ���������ݸ���
		int *size1 = new int[k];
		//������ĸ��º�ÿ���������ݸ���
		int *size2 = new int[k];
		//�Դ�����ĸ�����ʼ��
		for (i = 0; i < k; i++){
			size1[i] = 0;
		}
		//���ÿ�θ���ÿ�����ݱ仯ֵ
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
			//TODO���������
		}
		int count2 = 0;
		int count3 = 0;





		//���ķ���
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
			//���ϴ���ֹͬͣѭ��
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
		//���һ������ĸ���
		//int finalGroupCluster = Groupcount == 1 ? k : Groupcount % tan == 0 ? tan : k % tan;
		//���ÿ�����ݶ�Ӧÿ����½�
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

		//���ȫ���Ͻ� 
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
					//�����ǰ�����ڵ���ʱ�������j������Ϊͬһ��Ⱥ��˵����һ�μ������Сֵ��ʱΪ�ڶ�С����˵����һ��Ⱥ���Ѿ�����õ㲻��������һ��Ⱥ�齫��һ���μ������Сֵ��Ϊ��һ��Ⱥ����½� 
					if (listclusterIndexNew[listclusterIndex[i]] == listclusterIndexNew[j]){
						secondDis = dis;
					}
					else{
						listsecondminumgroup[i][listclusterIndexNew[listclusterIndex[i]]] = dis;
					}
					dis = d;
					listclusterIndex[i] = j;

				}
				//���d����dis��С��secondDis˵����ʱ��dΪ�ڶ�Сֵ 
				else if (d<secondDis){
					secondDis = d;
				}
				//�����ǰ���Ѿ��Ǹ�Ⱥ������һ����ѵڶ�С�ĵ���Ϊ������½������������㲻�����Ⱥ��Ӧ�ð���Сֵ��Ϊ�½��ʱ������Ϊ������������һ��ѭ��������������½��и��� 
				if ((j + 1) == t1){
					listsecondminumgroup[i][top++] = secondDis;
					if (top1 != Groupcount){
						t1 = t1 + groupSize[top1++];
					}
					secondDis = DBL_MAX;
				}
			}
			//listsecondminum_distance[i] = secondDis;
			//�����ݵ��Ͻ����Ϊ��Сֵ 
			listsmallest_distance[i] = dis;
			size1[listclusterIndex[i]]++;
		}
		//�����仯���� 
		double groupbiggest_shift = -1;
		//ÿ������仯���� 
		double *listgroupbiggest_shift = new double[Groupcount];
		//ÿ�����ĵ�ı仯���� 
		double *listshift = new double[k];

		int top = 0, top1 = 1;
		f1 = true;
		//�ô�ͳ������������ 
		getMeanTwo(listclusterIndex, datas, arr, iniPoints, num, colum2, k, size);
		t1 = groupSize[0];
		for (i = 0; i < k; i++){
			//�������ĵ�仯���� 
			double shift = getDistance(iniPoints[i], arr[i], colum2);
			count1++;
			listshift[i] = shift;
			//�ҵ���Ⱥ������仯���� 
			groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
			//������������һ�����ĵ���Ҫ�����������벢�Ұ����仯�����ʼ�� 
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
		//���ϴ���ֹͬͣѭ��
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
		//������
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
				//����ǰ�����Ͻ���ϸĵ��������ĵı仯����֮����Ϊ�µ����Ͻ� 
				listsmallest_distance[i] = listsmallest_distance[i] + listshift[listclusterIndex[i]];
				for (j = 0; j < Groupcount; j++){
					listsecondminumgroup_old[j] = listsecondminumgroup[i][j];
				}
				for (j = 0; j < Groupcount; j++){
					//����һ�ε����½������仯����֮����Ϊ�µ����½� 
					listsecondminumgroup[i][j] -= listgroupbiggest_shift[j];
				}
				int listclusterIndex_old = listclusterIndex[i];
				//int indd = listsmallest_distance[i];
				//���ȫ�����½��Ѿ��������Ͻ�˵���ĵ㲻��Ҫ�ı���������ѭ������Ļ���Ҫ�������� 
				if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
					//�����Ͻ���СΪ�ĵ㵽��ǰ���ĵľ���
					double ubTight = getDistance(datas[i], iniPoints[listclusterIndex[i]], colum2);
					count1++;
					listsmallest_distance[i] = ubTight;
					//���ȫ���½�����С����С����Ͻ���ôҪѡ��ÿһ��С���Ͻ�������ĵ��п��ܴ������� 
					if (min(listsecondminumgroup[i], Groupcount) < listsmallest_distance[i]){
						for (j = 0; j < Groupcount; j++){

							if (listsecondminumgroup[i][j] < (listsmallest_distance[i])){
								int qqq = groupSize[j];
								noLocalCount += (double)qqq;
								//�����½����Ϊ������ 
								listsecondminumgroup[i][j] = DBL_MAX;
								listg[top++] = j;
							}
						}
						flag = listclusterIndex[i];
						for (j = 0; j < top; j++){
							int t = listg[j];
							int A = groupSize[t];

							for (m = 0; m < A; m++){
								//���ڵ�ǰ��СֵΪ�õ㵱ǰ�ľ����ľ�������������뵱ǰ��Сֵ�Ƚ������������������и������Բ��ض������бȽ� 
								if (groupPointsData[t][m] != listclusterIndex_old){
									/*	dis = DBL_MAX;
									d = getDistance(iniPoints[tan * t + m], datas[i], colum2);
									if (d < dis){
									dis = d;
									listclusterIndex[i] = tan * t + m;
									}*/
									//��������� 
									if (listsecondminumgroup[i][t] >listsecondminumgroup_old[t] - listshift[groupPointsData[t][m]]){

										double lbtight = getDistance(datas[i], iniPoints[groupPointsData[t][m]], colum2);
										localCount++;
										if (lbtight < listsmallest_distance[i]){
											//�����û�б���ô����ڶ�СҪ�����һ�ε���Сֵ��Ϊ��ǰ��Ϊ�������ڸ����½���ҪΪ�ڶ�С 
											if (listclusterIndexNew[flag] == t){
												//if (listclusterIndex[i] / tan == (tan * t + m) / tan){
												listsecondminumgroup[i][t] = listsmallest_distance[i];
											}
											else{
												//listsecondminumgroup[i][listclusterIndex[i] / tan] = listsmallest_distance[i];
												//listsecondminumgroup[i][flag / tan] = listsmallest_distance[i];
												listsecondminumgroup[i][listclusterIndexNew[flag]] = listsmallest_distance[i];
											}
											//�����Ͻ� 
											listsmallest_distance[i] = lbtight;
											//���´������� 
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
						//�˴�Ҫ�Ա仯�����ݲ�ֵ����ͳ�����Ż��������ĵ��㷨 
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
				//�������ĵ�仯���� 
				double shift = getDistance(iniPoints[i], arr[i], colum2);
				count1++;
				listshift[i] = shift;
				//�ҵ���Ⱥ������仯���� 
				groupbiggest_shift = groupbiggest_shift < shift ? shift : groupbiggest_shift;
				//������������һ�����ĵ���Ҫ�����������벢�Ұ����仯�����ʼ�� 
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
			//���ϴ���ֹͬͣѭ��
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




