#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <boost/random.hpp>


#define maxgen 1000
#define sizepop 100 // ��Ⱥ��Ŀ
//#define pcross 0.7 // �������
//#define pmutation 0.06 // �������
#define lenchrom 34 // Ⱦɫ�峤��(���ＴΪ���и���)

//int sizepop;
double pcross, pmutation;

double city_pos[lenchrom][2];
int chrom[sizepop][lenchrom]; // ��Ⱥ
int best_result[lenchrom]; // ���·��
double min_distance; // ���·������
int flag;
FILE* txtFile;

void init(void); // ��Ⱥ��ʼ������
double distance(double*, double*); // ������������֮��ľ���
double* min(double*); // ��������������Сֵ
double path_len(int*); // ����ĳһ��������·�����ȣ���Ӧ�Ⱥ���Ϊ·�߳��ȵĵ���
void Choice(int[sizepop][lenchrom]); // ѡ�����
void Cross(int[sizepop][lenchrom]); // �������
void Mutation(int[sizepop][lenchrom]); // �������
void Reverse(int[sizepop][lenchrom]); // ��ת����
void Readdata();
void Writeresult();

//void Writeresult(int i)
//{
//	FILE* fp2 = fopen("result.txt", "w");
//	fprintf(fp2, "This is %d genes: \n", i);
//	fprintf(fp2, "This is the best way for now: ");
//	for (int k = 0; k < lenchrom - 1; k++)
//	{
//		fprintf(fp2, "%d-->", best_result[k]);
//	}
//	fprintf(fp2, "%d \n", best_result[lenchrom - 1]);
//}
void Readdata()
{
	int i, j;
	FILE* fp1 = fopen("data.txt", "r");
	for (i = 0;i < lenchrom;i++)
		for (j = 0; j < 2; j++)
		{
			fscanf(fp1, "%lf", &city_pos[i][j]);
		}
}

//(0,1)�������������  
double randx()
{
	double val;

	boost::mt19937 generator(time(0) * rand());
	boost::uniform_real<> uniform_real_generate_x(0, 1);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

	val = random_real_num_x();

	return(val);
}

// ��Ⱥ��ʼ��
void init(void)
{
	int num = 0;
	while (num < sizepop)
	{
		for (int i = 0; i < sizepop; i++)
			for (int j = 0; j < lenchrom; j++)
				chrom[i][j] = j + 1;
		num++;
		for (int i = 0; i < lenchrom - 1; i++)
		{
			for (int j = i + 1; j < lenchrom; j++)
			{
				int temp = chrom[num][i];
				chrom[num][i] = chrom[num][j];
				chrom[num][j] = temp; // ������num������ĵ�i��Ԫ�غ͵�j��Ԫ��
				num++;
				if (num >= sizepop)
					break;
			}
			if (num >= sizepop)
				break;
		}
		// ������������ѭ�������޷������㹻�ĳ�ʼ���壬������ٲ���һ����
		// ���巽ʽ����ѡ����������λ�ã�Ȼ�󽻻�
		while (num < sizepop)
		{
			double r1 = randx();
			double r2 = randx();
		
			int p1 = (int)(lenchrom * r1); // λ��1
			int p2 = (int)(lenchrom * r2); // λ��2
			int temp = chrom[num][p1];
			chrom[num][p1] = chrom[num][p2];
			chrom[num][p2] = temp;    // ��������λ��
			num++;
		}
	}
}

// ���뺯��
double distance(double* city1, double* city2)
{
	double x1 = *city1;
	double y1 = *(city1 + 1);
	double x2 = *(city2);
	double y2 = *(city2 + 1);
	double dis = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return dis;
}
// min()����
double* min(double* arr)
{
	static double best_index[2];
	double min_dis = *arr;
	double min_index = 0;
	for (int i = 1; i < sizepop; i++)
	{
		double dis = *(arr + i);
		if (dis < min_dis)
		{
			min_dis = dis;
			min_index = i;
		}
	}
	best_index[0] = min_index;
	best_index[1] = min_dis;
	return best_index;
}

// ����·������
double path_len(int* arr)
{
	double path = 0; // ��ʼ��·������
	int index = *arr; // ��λ����һ������(�������)
	for (int i = 0; i < lenchrom - 1; i++)
	{
		int index1 = *(arr + i);
		int index2 = *(arr + i + 1);
		double dis = distance(city_pos[index1 - 1], city_pos[index2 - 1]);
		path += dis;
	}
	int last_index = *(arr + lenchrom - 1); // ���һ���������
	int first_index = *arr; // ��һ���������
	double last_dis = distance(city_pos[last_index - 1], city_pos[first_index - 1]);
	path = path + last_dis;
	return path; // �����ܵ�·������
}

// ѡ�����
void Choice(int chrom[sizepop][lenchrom])
{
	double pick;
	double choice_arr[sizepop][lenchrom];
	double fit_pro[sizepop];
	double sum = 0;
	double fit[sizepop]; // ��Ӧ�Ⱥ�������(����ĵ���)
	for (int j = 0; j < sizepop; j++)
	{
		double path = path_len(chrom[j]);
		double fitness = 1 / path;
		fit[j] = fitness;
		sum += fitness;
	}
	for (int j = 0; j < sizepop; j++)
	{
		fit_pro[j] = fit[j] / sum; // ��������
	}
	// ��ʼ���̶�
	for (int i = 0; i < sizepop; i++)
	{
		pick = randx(); // 0��1֮��������
		for (int j = 0; j < sizepop; j++)
		{
			pick = pick - fit_pro[j];
			if (pick <= 0)
			{
				for (int k = 0; k < lenchrom; k++)
					choice_arr[i][k] = chrom[j][k]; // ѡ��һ������
				break;
			}
		}

	}
	for (int i = 0; i < sizepop; i++)
	{
		for (int j = 0; j < lenchrom; j++)
			chrom[i][j] = choice_arr[i][j];
	}
}

//�������
void Cross(int chrom[sizepop][lenchrom])
{
	double pick;
	double pick1, pick2;
	int choice1, choice2;
	int pos1, pos2;
	int temp;
	int conflict1[lenchrom]; // ��ͻλ��
	int conflict2[lenchrom];
	int num1, num2;
	int index1, index2;
	int move = 0; // ��ǰ�ƶ���λ��
	while (move < lenchrom - 1)
	{
		pick = randx(); // ���ھ����Ƿ���н������
		if (pick > pcross)
		{
			move += 2;
			continue; // ���β����н���
		}
		// ���ò���ӳ���ӽ�
		choice1 = move; // ����ѡȡ�ӽ�����������
		choice2 = move + 1; // ע������±�Խ��
		pick1 = randx();
		pick2 = randx();
		pos1 = (int)(pick1 * lenchrom); // ����ȷ�������ӽ����λ��
		pos2 = (int)(pick2 * lenchrom);
		while (pos1 > lenchrom - 2 || pos1 < 1)
		{
			pick1 = randx();
			pos1 = (int)(pick1 * lenchrom);
		}
		while (pos2 > lenchrom - 2 || pos2 < 1)
		{
			pick2 = randx();
			pos2 = (int)(pick2 * lenchrom);
		}
		if (pos1 > pos2)
		{
			temp = pos1;
			pos1 = pos2;
			pos2 = temp; // ����pos1��pos2��λ��
		}
		for (int j = pos1; j <= pos2; j++)
		{
			temp = chrom[choice1][j];
			chrom[choice1][j] = chrom[choice2][j];
			chrom[choice2][j] = temp; // �������˳��
		}
		num1 = 0;
		num2 = 0;
		if (pos1 > 0 && pos2 < lenchrom - 1)
		{
			for (int j = 0; j <= pos1 - 1; j++)
			{
				for (int k = pos1; k <= pos2; k++)
				{
					if (chrom[choice1][j] == chrom[choice1][k])
					{
						conflict1[num1] = j;
						num1++;
					}
					if (chrom[choice2][j] == chrom[choice2][k])
					{
						conflict2[num2] = j;
						num2++;
					}
				}
			}
			for (int j = pos2 + 1; j < lenchrom; j++)
			{
				for (int k = pos1; k <= pos2; k++)
				{
					if (chrom[choice1][j] == chrom[choice1][k])
					{
						conflict1[num1] = j;
						num1++;
					}
					if (chrom[choice2][j] == chrom[choice2][k])
					{
						conflict2[num2] = j;
						num2++;
					}
				}

			}
		}
		if ((num1 == num2) && num1 > 0)
		{
			for (int j = 0; j < num1; j++)
			{
				index1 = conflict1[j];
				index2 = conflict2[j];
				temp = chrom[choice1][index1]; // ������ͻ��λ��
				chrom[choice1][index1] = chrom[choice2][index2];
				chrom[choice2][index2] = temp;
			}
		}
		move += 2;
	}
}

// �������
// ������Բ�ȡ���ѡȡ�����㣬����Ի�λ��
void Mutation(int chrom[sizepop][lenchrom])
{
	double pick, pick1, pick2;
	int pos1, pos2, temp;
	for (int i = 0; i < sizepop; i++)
	{
		pick = randx(); // �����ж��Ƿ���б������
		if (pick > pmutation)
			continue;
		pick1 = randx();
		pick2 = randx();
		pos1 = (int)(pick1 * lenchrom); // ѡȡ���б����λ��
		pos2 = (int)(pick2 * lenchrom);
		while (pos1 > lenchrom - 1)
		{
			pick1 = randx();
			pos1 = (int)(pick1 * lenchrom);
		}
		while (pos2 > lenchrom - 1)
		{
			pick2 = randx();
			pos2 = (int)(pick2 * lenchrom);
		}
		temp = chrom[i][pos1];
		chrom[i][pos1] = chrom[i][pos2];
		chrom[i][pos2] = temp;
	}
}

void Reverse(int chrom[sizepop][lenchrom])
{
	double pick1, pick2;
	double dis, reverse_dis;
	int n;
	int flag, pos1, pos2, temp;
	int reverse_arr[lenchrom];

	for (int i = 0; i < sizepop; i++)
	{
		flag = 0; // ���ڿ��Ʊ�����ת�Ƿ���Ч
		while (flag == 0)
		{
			pick1 = randx();
			pick2 = randx();
			pos1 = (int)(pick1 * lenchrom); // ѡȡ������ת������λ��
			pos2 = (int)(pick2 * lenchrom);
			while (pos1 > lenchrom - 1)
			{
				pick1 = randx();
				pos1 = (int)(pick1 * lenchrom);
			}
			while (pos2 > lenchrom - 1)
			{
				pick2 = randx();
				pos2 = (int)(pick2 * lenchrom);
			}
			if (pos1 > pos2)
			{
				temp = pos1;
				pos1 = pos2;
				pos2 = temp; // ����ʹ��pos1 <= pos2
			}
			if (pos1 < pos2)
			{
				for (int j = 0; j < lenchrom; j++)
					reverse_arr[j] = chrom[i][j]; // ��������
				n = 0; // ��ת��Ŀ
				for (int j = pos1; j <= pos2; j++)
				{
					reverse_arr[j] = chrom[i][pos2 - n]; // ��ת����
					n++;
				}
				reverse_dis = path_len(reverse_arr); // ��ת֮��ľ���
				dis = path_len(chrom[i]); // ԭʼ����
				if (reverse_dis < dis)
				{
					for (int j = 0; j < lenchrom; j++)
						chrom[i][j] = reverse_arr[j]; // ���¸���
				}
			}
			flag = 1;
		}

	}
}

// ������
void start_function(void)
{
	char result_file[100];

	//���txt�ļ��������
	if (flag == 0)
	{
		sprintf(result_file, "result/result_p.txt");
			if ((txtFile = fopen(result_file, "a")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		else
		{
			if ((txtFile = fopen(result_file, "a")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		}
	}
	else if (flag == 1) {
		sprintf(result_file, "result/result_pc.txt");

		if (pcross == 0.1) {
			if ((txtFile = fopen(result_file, "w")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		}
		else
		{
			if ((txtFile = fopen(result_file, "a")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		}
	}
	else if (flag == 2)
	{
		sprintf(result_file, "result/result_pm.txt");

		if (pmutation == 0.01) {
			if ((txtFile = fopen(result_file, "w")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		}
		else
		{
			if ((txtFile = fopen(result_file, "a")) == NULL)
			{
				printf("�ļ���ʧ�ܣ�");
				exit(1);
			}
		}
	}
	else
	{
		printf("error\n");
		exit(1);
	}

	printf("\n F���ڼ��㣬���Ժ�...\n");

	//FILE* fp2 = fopen("result.txt", "w");
	Readdata();
	//int maxgen;
	//printf("Please input the genes:\n");
	//scanf("%d", &maxgen);
	time_t start, finish;
	start = clock(); // ��ʼ��ʱ
	srand((unsigned)time(NULL)); // ��ʼ�����������
	init(); // ��ʼ����Ⱥ

	int best_fit_index = 0; //���·�����ִ���
	double distance_arr[sizepop];
	double dis;
	for (int j = 0; j < sizepop; j++)
	{
		dis = path_len(chrom[j]);
		distance_arr[j] = dis;
	}
	double* best_index = min(distance_arr); // �������·�������
	min_distance = *(best_index + 1); // ���·��
	int index = (int)(*best_index); // ���·�����
	for (int j = 0; j < lenchrom; j++)
		best_result[j] = chrom[index][j]; // ���·������

										  // ��ʼ����
	double* new_arr;
	double new_min_dis;
	int new_index;
	for (int i = 0; i < maxgen; i++)
	{
		Choice(chrom); // ѡ��
		Cross(chrom); //����
		Mutation(chrom); //����
		Reverse(chrom); // ��ת����
		for (int j = 0; j < sizepop; j++)
			distance_arr[j] = path_len(chrom[j]); // ��������
		new_arr = min(distance_arr);
		new_min_dis = *(new_arr + 1); //�µ����·��
		if (new_min_dis < min_distance)
		{
			min_distance = new_min_dis; // �������·��
			new_index = (int)(*new_arr);
			for (int j = 0; j < lenchrom; j++)
				best_result[j] = chrom[new_index][j]; // �������·������
			best_fit_index = i + 1; // ���·������
		}

		if (flag == 0)
			fprintf(txtFile, "%d %d %.2f\n", sizepop, i+1, min_distance);
		else if (flag == 1)
			fprintf(txtFile, "%.2f %d %.2f\n", pcross, i+1, min_distance);
		else if (flag == 2)
			fprintf(txtFile, "%.2f %d %.2f\n", pmutation, i+1, min_distance);
		else
		{
			printf("error\n");
			exit(1);
		}
		//fprintf(txtFile, "%d %lf\n", i + 1, min_distance);

	//	fprintf(fp2, "This is %d genes: \n", i);
	//	fprintf(fp2, "This is the best way for now: ");
	//	for (int k = 0; k < lenchrom - 1; k++)
	//	{
	//		fprintf(fp2, "%d-->", best_result[k]);
	//	}
	//	fprintf(fp2, "%d \n", best_result[lenchrom - 1]);
	}
	//finish = clock(); // �������
	//double duration = ((double)(finish - start)) / CLOCKS_PER_SEC; // �����ʱ
	//fprintf(fp2, "This is solve %d cities' Tsp,population number:%d,evolution number:%d\n", lenchrom, sizepop, maxgen);
	//for (int k = 0; k < lenchrom - 1; k++)
	//{
	//	fprintf(fp2, "%d-->", best_result[k]);
	//}
	//fprintf(fp2, "%d \n", best_result[lenchrom - 1]);
	//fprintf(fp2, "The length of best way: %lf,get the result in %d genes\n", min_distance, best_fit_index);
	//fprintf(fp2, "Running time: %lf seconds.\n", duration);
	fclose(txtFile);
	/*getchar();
	system("pause");*/
	printf("����������result.txt\n");
	
}

//������Ⱥ��С���㷨������Ӱ��
//POPSIZE [20,200]
//PXOVER = 0.7 , PMUTATION = 0.06
void test_popsize(double pc, double pm)
{
	int p;

	pcross = pc;
	pmutation = pm;

	flag = 0; //�������
	for (p = 10; p <= 100; p += 10)
	{
		//sizepop = p;
		printf("\n\n p = %d, pm = %lf, pc = %lf\n", sizepop, pmutation, pcross);

		//���������������ֵ
	   /* for (int j = 0; j < POPSIZE; j++)
		{
			for (int i = 0; i < NVARS; i++)
			{
				printf("P = %d, N = %d, x%d = %lf\n", j+1,i+1, i+1,population[j].gene[i]);
			}
		}*/

		//�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
		start_function(); //Rotated High Conditioned Elliptic Funcktion


		//�ͷ��ڴ�ռ�
		/*free(x);
		free(f);
		delete[]population;
		delete[]newpopulation;*/
	}
}


//���Խ����С���㷨������Ӱ�� 
//PXOVER [0.9.0.1]
//POPSIZE = 50, PMUTATION = 0.06
void test_pxover(int p, double pm)
{
	double pc;

	//sizepop = p;
	pmutation = pm;

	flag = 1; //�������

	for (pc = 0.1; pc <= 1.0; pc += 0.1)
	{
		pcross = pc;
		printf("\n\n p = %d, pm = %lf, pc = %lf\n", sizepop, pmutation, pcross);

		//���������������ֵ
		/* for (int j = 0; j < POPSIZE; j++)
		 {
			 for (int i = 0; i < NVARS; i++)
			 {
				 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
			 }
		 }*/

		 //�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
		start_function(); //Rotated High Conditioned Elliptic Funcktion



	}

	//�ͷ��ڴ�ռ�
	//free(x);
	//free(f);
	//delete[]population;
	//delete[]newpopulation;
}


//���Ա����С���㷨������Ӱ��
//PMUTATION [0.01.0.1]
//POPSIZE = 50, PXOVER = 0.7
void test_pmutation(int p, double pc)
{
	double pm;

	//sizepop = p;
	pcross = pc;

	flag = 2; //�������

	for (pm = 0.01; pm <= 0.1; pm += 0.01)
	{
		pmutation = pm;
		printf("\n\n p = %d, pm = %lf, pc = %lf\n", sizepop, pmutation, pcross);


		//���������������ֵ
		/* for (int j = 0; j < POPSIZE; j++)
		 {
			 for (int i = 0; i < NVARS; i++)
			 {
				 printf("P = %d, N = %d, x%d = %lf\n", j + 1, i + 1, i + 1, population[j].gene[i]);
			 }
		 }*/

		 //�������Ժ������粻����ĳ���Ժ�����ͨ��"//"ע�͵�����
		start_function(); //Rotated High Conditioned Elliptic Funcktion

	}

	//�ͷ��ڴ�ռ�
	/*free(x);
	free(f);
	delete[]population;
	delete[]newpopulation;*/
}


//������
void main(void) {

	int p;
	double pc, pm;

	p = 100;   //��Ⱥ��С [20,200]
	pc = 0.7;   //������  [0.9,0.1]
	pm = 0.06;   //������ [0.01,0.1]

	//flag = 0;
	//pcross = 0.7;
	//pmutation = 0.06;
	//test_popsize(pc, pm); //������Ⱥ��С���㷨������Ӱ��

	test_pxover(p, pm); //���Խ����С���㷨������Ӱ��

	test_pmutation(p, pc); //���Ա����С���㷨������Ӱ��

	//start_function();

}