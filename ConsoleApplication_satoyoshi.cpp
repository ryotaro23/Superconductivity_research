/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
本プログラムについて

作成日：2017年7月9日
作成者：橋本 尚樹


本プログラムは,諸々の計算を別々のプログラムで書いてきたことによる,プログラムの散在化を
今以上に進行することを防ぐために作るものである.
ゆえに,本プログラムは研究に必要となる,あらゆる計算に対し,柔軟に耐えうるものでなければならない.


引継ぎ：2020年 大科智樹
引継ぎ：2021年 佐藤吉祥
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
更新内容記録
SENTENSE(XX年YY月ZZ日)
・出力ファイルに,無次元の数値に対応した,実際の物理量も一緒に出すようにした(2017年8月22日)
・OUTPUT_CD_VIEW_FILE を追加.より簡単に,ファイルを出力・不出力を選択できるようにするため(2017年8月23日)


・2ndMF、1st+1、1st-1について導入(2021.10.15)

・周期2のアンチドットについて対応

・アニール時間の設定導入(2021.12.15)

・ピニング力についてカットオフ設定範囲を長方形に(2022.02.11)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*プリプロセッサ　ヘッダファイル定義															*/	
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
//#include <Windows.h>
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_pow_int.h"
#include <direct.h>
#include <time.h>
//#include "rand.cpp"
#include "rand_1 2.cpp"




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*プリプロセッサ　定数定義																		*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#define PINNING_SITE //ピニングサイトを考慮するか
//#define Fp_INSIDE  //ピニングサイト内のポテンシャルを考慮するか
#define LORENTZ      //ローレンツ力を考慮するか
#define THERMAL      //熱による揺らぎの効果を考慮するか

//#define ANNEAL //アニールするか
#ifdef ANNEAL
#define ANNEAL_TIME 700000  //アニールをどの程度の時間(step数行うか)
#endif


#define OUTPUT_CD_VIEW_FILE //cd_viewファイルが要らなければ,コメントアウトしてください

#define VARIABLE     //変数を導入するか
#define CHANGE_KLX   //ローレンツ力を変数として動かすか
//#define CHANGE_DP  //ピニングサイト間距離を変数として動かすか
//#define CHANGE_KT  //温度を変数として動かすか
//#define CHANGE_TPM //ローレンツ力の周期を変数として動かすか

//#define SQUARE_CD   //ボルテックスの初期配置を四角格子にするか
//#define TRIANGLE_CD //ボルテックスの初期配置を三角格子にするか
//#define PS_EDGE_CD    //ボルテックスの初期配置をピニングサイトの左端に合わせるか(低ローレンツ領域での見かけ上のボルテックス平均速度の有限の値をゼロに合わせるため)
//#define RANDOM_CD   //ボルテックスの初期配置をランダムにするか

//#define PS_2ND_MF_CD
//#define PS_2ND_MF_CD_2PAIR

//#define PS_LAST_CD_RAND			//最後の1本のボルテックス配置をランダムにするか
//#define PS_LAST_CD_DISIGNATE		//最後の1本のボルテックス配置を指定するか

#define PS_DUO_CD


//#define SET_PINNING_SITE_3
//#define SET_PINNING_SITE_3_2PAIR
#define SET_PINNING_SITE_2


#define N		12//7//9//いくつのボルテックスを使うか
#define Np		12//9//いくつのピニングサイトを使うか
#define sqrNp	2//(int)pow(N,0.5)   //周期3のときは3?
#define sqrN	2
#define NX		3
#define NY		3

#define PP		(Np/3)
//#define STOP 1.37e-17	//struve_H1用
#define K1(x) gsl_sf_bessel_K1(x)	//1次の第二種ベッセル


#define POS_OFFS 0

//char PROGRAM_NAME[64];
//sprintf(PROGRAM_NAME,"VRE_MD");
#define PROGRAM_NAME "VRE_MD"	//プログラム名
#define PROGRAM_NUM	140	//プログラムの番号
//#define save_num	50		//ファイルナンバー

#define PC "DELL_INSPIRON_15"	//プログラムがあるPCの名前
#define SAVE_DATA "DOCUMENT"	//プログラムがあるフォルダの名前




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*プロトタイプ宣言																				*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void output_file(double* cd,double* ps, double*fc ,double*fc_vvi_p,double*fc_vvi_m ,double*fc_pin,double*fc_Lor,double*fc_tho,double side, double side_x, double side_y, int save_num);
void calc_f(double* cd, double* fc, double*fc_vvi_p,double*fc_vvi_m,double*fc_pin,double*fc_Lor,double*fc_thm,double* vl, double* ps, double side, double side_x, double side_y);
void calc_potential(double* ps, double side);

void set_pinning_site_3(double* ps, double side, double Rl, double Rm, double Rs);
void set_pinning_site_3_2Pair(double* ps, double side, double Rl, double Rm, double Rs);
void set_pinning_site_2(double* ps, double side, double Rl, double Rs);


double fvv(const double r);
double fp(const double r);
double fp_inside(double r);
double struve_H1(const double x);
double mean(double data[], int step, int n);
double f(double* vl, double* fc, int j);

void set_cd_square_and_vl(double* cd, double* vl, double side, double side_x, double side_y);
void set_cd_triangle_and_vl(double* cd, double* vl, double side, double side_x, double side_y);
void set_cd_psedge_and_vl(double* cd, double* vl, double side, double side_x, double side_y);
void set_cd_random_and_vl(double* cd, double* vl, double side, double side_x, double side_y);

//2021.10.22、セカンドマッチングフィールド考慮
void set_cd_ps_2ndMF_and_vl(double* cd, double* vl, double side, double side_x, double side_y,double rl,double rm,double rs);
void set_cd_ps_2ndMF_2Pair_and_vl(double* cd, double* vl, double side, double side_x, double side_y, double rl, double rm, double rs);

//2021.11.11  周期2のピニングサイトについて考慮
void set_cd_ps_duo_and_vl(double* cd, double* vl, double side, double side_x, double side_y, double rl, double rm, double rs);


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*変数宣言・定義																				*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double
d = 1.0,			//膜厚
f0 = 1.0,			//VV相互作用の大きさの係数
Eta = 1.0,			//粘性係数
Kp = 2.0,			//ピン止め力係数
//Lp=0.05,
Lp = 0.3 * sqrt(2.0),	//ピニング回復長
Klx = 1.965,			//x方向ローレンツ力係数
Kly = 0.0,			//y方向ローレンツ力係数
lambda = 1.0,
//lambda = 10.0,
Lambda = 2 * lambda / tanh(d / lambda),
//Dp = 4.0,
Dp = 16,		// 実際の系だとDp=40(4um)

//2021追加分
/*PSLargeとPSMiddle間の中心間距離。PSMiddleとPSSmall間の中心間距離。PSSmallとPSLarge間の中心間距離。*/
Dp_L2M = 10,
Dp_M2S = 10,
Dp_S2L = 10,


//Dp=20,
Rl = 1.5,           //Large PSの半径
Rm = 1.0,
Rs = 0.5,
Fp_cutoff = Dp;//Dp;

double
fq = 2.068e-15,
mu_vac = 4 * M_PI * pow(10.0, -7.0),
boltzmann = 1.381e-23;

double    variable;
//double    set_temp = 8.4;      //実際の温度を何[K]に設定するか
double    set_temp = 4.2;

double	  dt = 0.001;
//int       period = 10000000;
//int       period = 100000;      //ローレンツ力の周期 20Mhz
int		  period = 25000;		//5Mhz
int       change_pm = period / 2; //ローレンツ力の±を切り替える 50000 (periodが半分まで到達したら±切り替える)
int       N_period = 3;        //交流で何周期分計算するか(左右にローレンツ力を何回振るか)
int       total_step = N_period * period;

//const int PP = 4;//ピニングサイトの周期


double
scale_length = 100e-9,
real_lambda = scale_length * lambda,
real_d = scale_length * d,
scale_visco = 2.66e-9,//この値は,結晶のNbを想定している可能性がある
scale_force = pow(fq, 2.0) / (2 * M_PI * mu_vac * pow(real_lambda, 3.0)),
scale_time = scale_length * scale_visco / scale_force,
scale_velocity = scale_length / scale_time,
scale_temp = pow(fq, 2.0) * real_d * dt / (4 * M_PI * mu_vac * boltzmann * pow(real_lambda, 2.0));

double
Kt = sqrt(4 * M_PI * mu_vac * boltzmann * pow(real_lambda, 2.0) * set_temp / (pow(fq, 2.0) * real_d * dt));

double
real_Lp = scale_length * Lp,
real_Dp = scale_length * Dp,
real_Rl = scale_length * Rl,
//real_Rm        =scale_length*Rm,
real_Rs = scale_length * Rs,
real_Fp_cutoff = scale_length * Fp_cutoff,

real_f0 = scale_force * f0,
Lorentz_x = scale_force * f0 * Klx,//x軸方向の単位長さあたりのローレンツ力
Lorentz_y = scale_force * f0 * Kly,

real_Eta = scale_visco * Eta,//この値は,結晶のNbを想定している可能性がある

real_total_time = scale_time * total_step * dt,
real_time_step = scale_time * dt,

temp = scale_temp * pow(Kt, 2.0);

rand_gauss* gauss = new rand_gauss();//

double  R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio;

double PinningMeshPos_X[Np][2];
double PinningMeshPos_Y[Np][2];
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*main関数																						*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int main()
{
	clock_t start, end;

	FILE* fpin;
	FILE* f_number_save;
	char fname[] = "variable_R_ratio.txt";

	//double  R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio;
	double Dp_S2L_min, Dp_S2L_max, dDp_S2L;

	int save_num;

	errno_t err;//errno_t型(int型)
	errno_t err_2;


	/****************************************************************************************/
	/*概要：variable_R_ratioファイル読み込み　												*/
	/*																						*/
	/*目的：パラメータ読み込み																*/
	/*																						*/
	/****************************************************************************************/
	err = fopen_s(&fpin, "variable_R_ratio.txt", "r");//fileを開く。失敗するとエラー
	if (err != 0) {
		printf("%s file not open!\n", "variable_R_ratio.txt");
		return err;
	}
	else {
		printf("%s file opened!\n", "variable_R_ratio.txt");
	}
	fscanf_s(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &R_standard_min, &R_standard_max, &dR_standard, &Rl_ratio, &Rm_ratio, &Rs_ratio,&Dp_L2M,&Dp_M2S,
														&Dp_S2L_min,&Dp_S2L_max,&dDp_S2L);
	//fscanf_s(fpin, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &R_standard_min, &R_standard_max, &dR_standard, &Rl_ratio, &Rm_ratio, &Rs_ratio, &Dp_L2M, &Dp_M2S, &Dp_S2L);
	//printf("%.3f %.3f %.3f %.3f %.3f %.3f\n", R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio);
    fclose(fpin);


	/****************************************************************************************/
	/*概要：num_data_saveファイル読み込み　ファイルナンバーを自動で記録可能					*/
	/*																						*/											
	/*目的：ファイル書き出しのための番号振り分けを自動で行うため。							*/
	/*		これを行わないと、前の計算ファイルに上書き(or追記)されてしまう					*/
	/*																						*/
	/****************************************************************************************/
	err_2 = fopen_s(&f_number_save, "num_data_save.txt", "r");//fileを開く。失敗するとエラー
	if (err != 0) {
		printf("%s file not open!\n", "num_data_save.txt");
		return err;
	}
	else {
		printf("%s file opened for reading!\n", "num_data_save.txt");
	}
	fscanf_s(f_number_save, "%d", &save_num);
	fclose(f_number_save);

	err_2 = fopen_s(&f_number_save, "num_data_save.txt", "w");//fileを開く。失敗するとエラー
	if (err != 0) {
		printf("%s file not open!\n", "num_data_save.txt");
		return err;
	}
	else {
		printf("%s file opened for writing!\n", "num_data_save.txt");
	}
	fprintf(f_number_save, "%d", save_num + 1);
	fclose(f_number_save);




	double real_Rm = scale_length * Rm;

#ifdef VARIABLE
	double R_standard=R_standard_min;
	//for (R_standard = R_standard_min; R_standard <= R_standard_max; R_standard += dR_standard)
	for (Dp_S2L = Dp_S2L_min; Dp_S2L <= Dp_S2L_max; Dp_S2L += dDp_S2L)
	{
		double
			Rl = R_standard * Rl_ratio,
			Rm = R_standard * Rm_ratio,
			Rs = R_standard * Rs_ratio,

			real_Rl = scale_length * Rl,
			real_Rm = scale_length * Rm,
			real_Rs = scale_length * Rs;

		printf("%.3f %.3f %.3f\n", Rl, Rm, Rs);
		printf("%.3f %.3f %.3f\n", R_standard_min, R_standard_max, dR_standard);
		printf("%.3f %.3f %.3f\n", Dp_L2M, Dp_M2S, Dp_S2L);

		//double min_variable = 1.000, max_variable = 2.070 + 0.0001, d_variable = 0.01;        //変数の定義域や刻み幅を設定 0.001
		double min_variable = 1.860, max_variable = 2.070 + 0.0001, d_variable = 0.001;
		//double min_variable = 1.780, max_variable = 1.850 + 0.0001, d_variable = 0.01;
		//double min_variable = 1.900, max_variable = 2.070 + 0.0001, d_variable = 0.001;        //変数の定義域や刻み幅を設定 0.001
#ifdef OUTPUT_CD_VIEW_FILE
		min_variable = 1.965;
		max_variable = 1.965 + 0.0001;
		//min_variable = 1.98;
		//max_variable = 1.98 + 0.0001;
		d_variable = 0.01;			//アニメーション用
#endif
		//d_variable = 0.01;
		for (variable = min_variable; variable <= max_variable; variable += d_variable)
		{
			start = clock();
#endif
			int i, j;
			int step;
			double cd[N * 3], vl[N * 3], fc[N * 3], fc0[N * 3], ps[Np * 3];
			double fc_vvi_p[N * 3], fc_vvi_m[N * 3];
			double fc_pin[N * 3];
			double fc_Lor[N * 3];
			double fc_tho[N * 3];
#ifdef PS_EDGE_CD
			double side = Dp * sqrN;
			//double side_x = Dp * NX;		//　計算で使用する一周期範囲X
			double side_x = Dp_L2M + Dp_M2S+Dp_S2L;
			double side_y = Dp * NY;		//　計算で使用する一周期範囲Y
			//double side_y = Dp_L2M + Dp_M2S + Dp_S2L;
			double sideh = side / 2;
			double sideh_x = side_x / 2.0;
			double sideh_y = side_y / 2.0;
#endif
#ifdef RANDOM_CD
			
			double side = Dp * sqrN;
			//double side_x = Dp * NX;		//　計算で使用する一周期範囲X
			//double side_x = Dp_L2M + Dp_M2S + Dp_S2L;
			double side_x = (Np / 6) * (Dp_L2M + Dp_M2S + Dp_S2L); //444でピニングサイト構成する際はこっち
			double side_y = Dp * NY;		//　計算で使用する一周期範囲Y
			//double side_y = Dp_L2M + Dp_M2S + Dp_S2L;
			double sideh = side / 2;
			double sideh_x = side_x / 2.0;
			double sideh_y = side_y / 2.0;
#endif
#ifdef SQUARE_CD
			double side = Dp * sqrN, side_x = Dp * NX, side_y = Dp * NY, sideh = side / 2, sideh_x = side_x / 2.0, sideh_y = side_y / 2.0;
#endif
#ifdef TRIANGLE_CD
			double side = Dp * sqrN, side_x = Dp * NX, side_y = Dp * NY, sideh = side / 2, sideh_x = side_x / 2.0, sideh_y = side_y / 2.0;
#endif
#ifdef PS_2ND_MF_CD

			double side = Dp * sqrN;
			//double side_x = Dp * NX;		//　計算で使用する一周期範囲X
			double side_x = Dp_L2M + Dp_M2S + Dp_S2L;
			double side_y = Dp * NY;		//　計算で使用する一周期範囲Y
			//double side_y = Dp_L2M + Dp_M2S + Dp_S2L;
			double sideh = side / 2;
			double sideh_x = side_x / 2.0;
			double sideh_y = side_y / 2.0;

#endif
#ifdef PS_2ND_MF_CD_2PAIR

			double side = Dp * sqrN;
			//double side_x = Dp * NX;		//　計算で使用する一周期範囲X
			double side_x = 2*(Dp_L2M + Dp_M2S + Dp_S2L);
			double side_y = Dp * NY;		//　計算で使用する一周期範囲Y
			//double side_y = Dp_L2M + Dp_M2S + Dp_S2L;
			double sideh = side / 2;
			double sideh_x = side_x / 2.0;
			double sideh_y = side_y / 2.0;

#endif
#ifdef PS_DUO_CD
			double side = Dp * sqrN;
			//double side_x = Dp * NX;		//　計算で使用する一周期範囲X
			//double side_x = Dp_L2M + Dp_M2S + Dp_S2L;
			double side_x = (Np/6)*(Dp_L2M + Dp_M2S + Dp_S2L); //444でピニングサイト構成する際はこっち
			double side_y = Dp * NY;		//　計算で使用する一周期範囲Y
			//double side_y = Dp_L2M + Dp_M2S + Dp_S2L;
			double sideh = side / 2;
			double sideh_x = side_x / 2.0;
			double sideh_y = side_y / 2.0;
#endif
			double cut_off2 = 36.0;
			double
				real_boxsize_x = scale_length * side_x,
				real_boxsize_y = scale_length * side_y;
			double V = 0.0;
			double V_plus = 0.0, V_minus = 0.0;//vlの＋とマイナスをそれぞれ計算

			double V_Row[Np / PP];		//各行ごとの速さ
			double V_Row_Plus[Np / PP], V_Row_Minus[Np / PP];

			double V_Col[PP];			//各列ごとの速さ
			double V_Col_Plus[PP], V_Col_Minus[PP];
#ifdef OUTPUT_CD_VIEW_FILE
			FILE* fp;
			char filename[128];
#ifdef VARIABLE
			char dir_1[64];
			char dir_2[64];
			char dir_3[64];
			sprintf_s(filename, "./%s_%d/%s_file_%d/variable_%f/%s_%d_%d_condition.txt"
				, PROGRAM_NAME, PROGRAM_NUM
				, PROGRAM_NAME, save_num
				, abs(variable)
				, PROGRAM_NAME, PROGRAM_NUM, save_num);
			sprintf_s(dir_1, "./%s_%d", PROGRAM_NAME, PROGRAM_NUM);
			sprintf_s(dir_2, "./%s_%d/%s_file_%d", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num);
			sprintf_s(dir_3, "./%s_%d/%s_file_%d/variable_%f", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, abs(variable));
			if (_mkdir(dir_1) == 0) { printf("SUCCESS\n"); }
			else { printf("FAILED\n"); }
			if (_mkdir(dir_2) == 0) { printf("SUCCESS\n"); }
			else { printf("FAILED\n"); }
			if (_mkdir(dir_3) == 0) { printf("SUCCESS\n"); }
			else { printf("FAILED\n"); }
#endif
#ifndef VARIABLE
			char dir_1[64];
			char dir_2[64];
			char dir_3[64];
			sprintf(filename, "./%s_%d/%s_file_%d/%s_%d_%d_condition.txt"
				, PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
			sprintf(dir_1, "./%s_%d", PROGRAM_NAME, PROGRAM_NUM);
			sprintf(dir_2, "./%s_%d/%s_file_%d", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num);
			if (_mkdir(dir_1) == 0) { printf("SUCCESS\n"); }
			else { printf("FAILED\n"); }
			if (_mkdir(dir_2) == 0) { printf("SUCCESS\n"); }
			else { printf("FAILED\n"); }
#endif
			errno_t err;
			if ((err = fopen_s(&fp, filename, "a")) == NULL)
			{
				fprintf(stderr, "file open error %s\n", fp);
			}
			fprintf(fp, "#PROGRAM_NAME :%s\n", PROGRAM_NAME);
			fprintf(fp, "#PROGRAM_NUM  :%d\n", PROGRAM_NUM);
			fprintf(fp, "#save_num     :%d\n", save_num);
			fprintf(fp, "#PC           :@%s\n", PC);
			fprintf(fp, "#SAVE_DATA    :@%s\n\n", SAVE_DATA);

			fprintf(fp, "#*#*#*#*#*#*#*シミュレーション条件#*#*#*#*#*#*#*\n");
#ifdef PINNING_SITE
			fprintf(fp, "#ピニングサイト                       :あり\n");
#endif
#ifndef PINNING_SITE
			fprintf(fp, "#ピニングサイト                       :なし\n");
#endif

#ifdef Fp_INSIDE
			fprintf(fp, "#ピニングサイト内でのポテンシャル勾配 :あり\n");
#endif
#ifndef Fp_INSIDE
			fprintf(fp, "#ピニングサイト内でのポテンシャル勾配 :なし\n");
#endif

#ifdef LORENTZ
			fprintf(fp, "#ローレンツ力                         :あり\n");
#endif
#ifndef LORENTZ
			fprintf(fp, "#ローレンツ力                         :なし\n");
#endif

#ifdef THERMAL
			fprintf(fp, "#熱散乱                               :あり\n");
#endif
#ifndef THERMAL
			fprintf(fp, "#熱散乱                               :なし\n");
#endif

#ifdef SQUARE_CD
			fprintf(fp, "#初期配置：正方格子配列\n");
#endif

#ifdef TRIANGLE_CD
			fprintf(fp, "#初期配置：正三角格子配列\n");
#endif

#ifdef PS_EDGE_CD
			fprintf(fp, "#初期配置：ピニングサイト左端に配置\n");
#endif

#ifdef RANDOM_CD
			fprintf(fp, "#初期配置：ランダム配置\n");
#endif

#ifdef VARIABLE

#ifdef CHANGE_KLX
			Klx = variable;
			Lorentz_x = scale_force * f0 * Klx;
			fprintf(fp, "#variable：Klx\n");
#endif 

#ifdef CHANGE_DP
			Dp = variable;
			fprintf(fp, "#variable：Dp\n");
#endif 

#ifdef CHANGE_KT
			set_temp = variable;
			Kt = sqrt(4 * M_PI * mu_vac * boltzmann * pow(real_lambda, 2.0) * set_temp / (pow(fq, 2.0) * real_d * dt));
			temp = scale_temp * pow(Kt, 2.0);
			fprintf(fp, "#variable：set_temp\n");
#endif 

#ifdef CHANGE_TPM
			change_pm = variable;
			fprintf(fp, "#variable：change_pm\n");
#endif 

#endif

#ifndef VARIABLE
			fprintf(fp, "#VARIABLEなし\n");
#endif

			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*#*長さ#*#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#膜厚                     d         %.3f | %.3e[m]\n", d, real_d);
			fprintf(fp, "#ピニング回復長           Lp        %.3f | %.3e[m]\n", Lp, real_Lp);
			fprintf(fp, "#磁場侵入長               lambda    %.3f | %.3e[m]\n", lambda, real_lambda);
			fprintf(fp, "#ピニングサイト中心間距離 Dp        %.3f | %.3e[m]\n", Dp, real_Dp);
			fprintf(fp, "#ピニングサイトL半径      Rl        %.3f | %.3e[m]\n", Rl, real_Rl);
			fprintf(fp, "#ピニングサイトM半径      Rm        %.3f | %.3e[m]\n", Rm, real_Rm);
			fprintf(fp, "#ピニングサイトS半径      Rs        %.3f | %.3e[m]\n", Rs, real_Rs);
			fprintf(fp, "#ピニングサイトS半径      Rs        %.3f | %.3e[m]\n", Rs, real_Rs);
			fprintf(fp, "#ピニング力カットオフ距離 Fp_cutoff %.3f | %.3e[m]\n", Fp_cutoff, real_Fp_cutoff);
			fprintf(fp, "#ボックスの大きさx軸方向  Lx        %.3f | %.3e[m]\n", side_x, real_boxsize_x);
			fprintf(fp, "#ボックスの大きさy軸方向  Ly        %.3f | %.3e[m]\n", side_y, real_boxsize_y);
			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*力の係数#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#力の単位                 f0  %.4f | %.4e[N/m]\n", f0, real_f0);
			fprintf(fp, "#粘性係数                 Eta %.4f | %.4e[kg/(m・s)]\n", Eta, real_Eta);
			fprintf(fp, "#ローレンツ力x方向        Klx %.4f | %.4e[N/m]\n", abs(Klx), abs(Lorentz_x));
			fprintf(fp, "#ローレンツ力y方向        Kly %.4f | %.4e[N/m]\n", abs(Kly), abs(Lorentz_y));
			fprintf(fp, "#温度                     Kt  %.4f | %.4e[K]\n", Kt, temp);
			fprintf(fp, "#ピニング力               Kp  %.4f | %.4e[N/m]\n", Kp, 2 * real_f0);
			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*その他#*#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#刻み時間      dt         %.1e | %.5e[sec]\n", dt, real_time_step);
			fprintf(fp, "#総ステップ数  total_step %d | %.5e[sec]\n", total_step, real_total_time);
#ifdef LORENTZ
			fprintf(fp, "#交流の周期    period     %d   | %.5e[sec]\n", period, scale_time * period);
			fprintf(fp, "#交流の周波数  frequency  %.5e[Hz]\n", 1 / (scale_time * period));
#endif
			fprintf(fp, "#VVIカットオフ cut_off2   %.3f | %.5e[m]\n", cut_off2, scale_length * cut_off2);
			fprintf(fp, "#粒子数        N          %d\n", N);
#ifdef VARIABLE
			fprintf(fp, "#variable下限 min_variable %f\n", min_variable);
			fprintf(fp, "#variable上限 max_variable %f\n", max_variable);
			fprintf(fp, "#variable刻幅 d_variable   %f\n", d_variable);
#endif
			fclose(fp);
#endif

#ifndef OUTPUT_CD_VIEW_FILE
			if (variable == min_variable || variable == max_variable)
				//if (variable >= min_variable && variable <= max_variable)		こちらの条件にすれば全パラメータでリアルスケール出力可能
			{
				FILE* fp;
				errno_t err;
				char filename[128];

#ifdef VARIABLE
				char dir_1[64];
				char dir_2[64];
				char dir_3[64];
				sprintf_s(filename, "./%s_%d/%s_file_%d/variable_%f/%s_%d_%d_condition.txt"
					, PROGRAM_NAME, PROGRAM_NUM
					, PROGRAM_NAME, save_num
					, fabs(variable)
					, PROGRAM_NAME, PROGRAM_NUM, save_num);
				sprintf_s(dir_1, "./%s_%d", PROGRAM_NAME, PROGRAM_NUM);
				sprintf_s(dir_2, "./%s_%d/%s_file_%d", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num);
				sprintf_s(dir_3, "./%s_%d/%s_file_%d/variable_%f", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, fabs(variable));
				if (_mkdir(dir_1) == 0) { printf("SUCCESS\n"); }
				else { printf("FAILED\n"); }
				if (_mkdir(dir_2) == 0) { printf("SUCCESS\n"); }
				else { printf("FAILED\n"); }
				if (_mkdir(dir_3) == 0) { printf("SUCCESS\n"); }
				else { printf("FAILED\n"); }
#endif
#ifndef VARIABLE
				char dir_1[64];
				char dir_2[64];
				char dir_3[64];
				sprintf(filename, "./%s_%d/%s_file_%d/%s_%d_%d_condition.txt"
					, PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
				sprintf(dir_1, "./%s_%d", PROGRAM_NAME, PROGRAM_NUM);
				sprintf(dir_2, "./%s_%d/%s_file_%d", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num);
				if (_mkdir(dir_1) == 0) { printf("SUCCESS\n"); }
				else { printf("FAILED\n"); }
				if (_mkdir(dir_2) == 0) { printf("SUCCESS\n"); }
				else { printf("FAILED\n"); }
#endif
				if ((err = fopen_s(&fp, filename, "wt")) == NULL)
				{
					fprintf(stderr, "file open error %s\n", filename);
				}

				fprintf(fp, "#PROGRAM_NAME :%s\n", PROGRAM_NAME);
				fprintf(fp, "#PROGRAM_NUM  :%d\n", PROGRAM_NUM);
				fprintf(fp, "#save_num     :%d\n", save_num);
				fprintf(fp, "#PC           :@%s\n", PC);
				fprintf(fp, "#SAVE_DATA    :@%s\n\n", SAVE_DATA);

				fprintf(fp, "#*#*#*#*#*#*#*シミュレーション条件#*#*#*#*#*#*#*\n");
#ifdef PINNING_SITE
				fprintf(fp, "#ピニングサイト                       :あり\n");
#endif
#ifndef PINNING_SITE
				fprintf(fp, "#ピニングサイト                       :なし\n");
#endif

#ifdef Fp_INSIDE
				fprintf(fp, "#ピニングサイト内でのポテンシャル勾配 :あり\n");
#endif
#ifndef Fp_INSIDE
				fprintf(fp, "#ピニングサイト内でのポテンシャル勾配 :なし\n");
#endif

#ifdef LORENTZ
				fprintf(fp, "#ローレンツ力                         :あり\n");
#endif
#ifndef LORENTZ
				fprintf(fp, "#ローレンツ力                         :なし\n");
#endif

#ifdef THERMAL
				fprintf(fp, "#熱散乱                               :あり\n");
#endif
#ifndef THERMAL
				fprintf(fp, "#熱散乱                               :なし\n");
#endif

#ifdef SQUARE_CD
				fprintf(fp, "#初期配置：正方格子配列\n");
#endif

#ifdef TRIANGLE_CD
				fprintf(fp, "#初期配置：正三角格子配列\n");
#endif

#ifdef PS_EDGE_CD
				fprintf(fp, "#初期配置：ピニングサイト左端に配置\n");
#endif

#ifdef RANDOM_CD
				fprintf(fp, "#初期配置：ランダム配置\n");
#endif

#ifdef VARIABLE

#ifdef CHANGE_KLX
				Klx = variable;
				Lorentz_x = scale_force * f0 * Klx;
				fprintf(fp, "#variable：Klx\n");
#endif 

#ifdef CHANGE_DP
				Dp = variable;
				fprintf(fp, "#variable：Dp\n");
#endif 

#ifdef CHANGE_KT
				set_temp = variable;
				Kt = sqrt(4 * M_PI * mu_vac * boltzmann * pow(real_lambda, 2.0) * set_temp / (pow(fq, 2.0) * real_d * dt));
				temp = scale_temp * pow(Kt, 2.0);
				fprintf(fp, "#variable：set_temp\n");
#endif 

#ifdef CHANGE_TPM
				change_pm = variable;
				fprintf(fp, "#variable：change_pm\n");
#endif 

#endif

#ifndef VARIABLE
				fprintf(fp, "#VARIABLEなし\n");
			}
#endif

			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*#*長さ#*#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#膜厚                     d         %.3f | %.3e[m]\n", d, real_d);
			fprintf(fp, "#ピニング回復長           Lp        %.3f | %.3e[m]\n", Lp, real_Lp);
			fprintf(fp, "#磁場侵入長               lambda    %.3f | %.3e[m]\n", lambda, real_lambda);
			fprintf(fp, "#ピニングサイト中心間距離 Dp        %.3f | %.3e[m]\n", Dp, real_Dp);
			fprintf(fp, "#ピニングサイトL半径      Rl        %.3f | %.3e[m]\n", Rl, real_Rl);
			fprintf(fp, "#ピニングサイトM半径      Rm        %.3f | %.3e[m]\n", Rm, real_Rm);
			fprintf(fp, "#ピニングサイトS半径      Rs        %.3f | %.3e[m]\n", Rs, real_Rs); 
			fprintf(fp, "#L~Mの距離                Dp_L2M    %.3f | %.3e[m]\n", Dp_L2M, Dp_L2M);  //6.30変更
			fprintf(fp, "#ピニング力カットオフ距離 Fp_cutoff %.3f | %.3e[m]\n", Fp_cutoff, real_Fp_cutoff);
			fprintf(fp, "#ボックスの大きさx軸方向  Lx        %.3f | %.3e[m]\n", side_x, real_boxsize_x);
			fprintf(fp, "#ボックスの大きさy軸方向  Ly        %.3f | %.3e[m]\n", side_y, real_boxsize_y);
			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*力の係数#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#力の単位                 f0  %.4f | %.4e[N/m]\n", f0, real_f0);
			fprintf(fp, "#粘性係数                 Eta %.4f | %.4e[kg/(m・s)]\n", Eta, real_Eta);
			fprintf(fp, "#ローレンツ力x方向        Klx %.4f | %.4e[N/m]\n", fabs(Klx), fabs(Lorentz_x));
			fprintf(fp, "#ローレンツ力y方向        Kly %.4f | %.4e[N/m]\n", fabs(Kly), fabs(Lorentz_y));
			fprintf(fp, "#温度                     Kt  %.4f | %.4e[K]\n", Kt, temp);
			fprintf(fp, "#ピニング力               Kp  %.4f | %.4e[N/m]\n", Kp, 2 * real_f0);
			fprintf(fp, "\n");

			fprintf(fp, "#*#*#*#*#*#*#*#*#*#*その他#*#*#*#*#*#*#*#*#*#*#*\n");
			fprintf(fp, "#刻み時間      dt         %.1e | %.5e[sec]\n", dt, real_time_step);
			fprintf(fp, "#総ステップ数  total_step %d | %.5e[sec]\n", total_step, real_total_time);
#ifdef LORENTZ
			fprintf(fp, "#交流の周期    period     %d   | %.5e[sec]\n", period, scale_time * period);
			fprintf(fp, "#交流の周波数  frequency  %.5e[Hz]\n", 1 / (scale_time * period));
#endif
			fprintf(fp, "#VVIカットオフ cut_off2   %.3f | %.5e[m]\n", cut_off2, scale_length * cut_off2);
			fprintf(fp, "#粒子数        N          %d\n", N);
#ifdef VARIABLE
			fprintf(fp, "#variable下限 min_variable %f\n", min_variable);
			fprintf(fp, "#variable上限 max_variable %f\n", max_variable);
			fprintf(fp, "#variable刻幅 d_variable   %f\n", d_variable);
#endif
			fclose(fp);
		}
#endif
		/****************************************************************************************/
		/*概要:ボルテックスの初期位置・初速度の設定												*/
		/*																						*/
		/****************************************************************************************/
		
#ifdef SQUARE_CD
		set_cd_square_and_vl(cd, vl, side, side_x, side_y);		//初期位置：四角格子
#endif
#ifdef TRIANGLE_CD
		set_cd_triangle_and_vl(cd, vl, side, side_x, side_y);	//初期位置：三角格子
#endif
#ifdef PS_EDGE_CD
		set_cd_psedge_and_vl(cd, vl, side, side_x, side_y);		//初期位置：ピニングサイト(PS)左端
#endif
#ifdef RANDOM_CD
		set_cd_random_and_vl(cd, vl, side, side_x, side_y);		//初期位置：ランダム  
#endif
#ifdef PS_2ND_MF_CD
		set_cd_ps_2ndMF_and_vl(cd, vl, side, side_x, side_y,Rl,Rm,Rs);   //初期位置：ピニングサイト(PS)指定位置。2ndMF考慮
#endif
#ifdef PS_2ND_MF_CD_2PAIR
		set_cd_ps_2ndMF_2Pair_and_vl(cd, vl, side, side_x, side_y, Rl, Rm, Rs);   //初期位置：ピニングサイト(PS)指定位置。2ndMF考慮
#endif
#ifdef PS_DUO_CD
		set_cd_ps_duo_and_vl(cd ,vl, side, side_x, side_y, Rl, Rm, Rs);
#endif
		/****************************************************************************************/
		/*概要:ボルテックス速さの配列変数初期化													*/
		/*																						*/
		/****************************************************************************************/

		for (i = 0; i < (Np/PP); i++)
		{
			V_Row[i] = 0;
			V_Row_Plus[i] = 0;
			V_Row_Minus[i] = 0;
		}
		for (i = 0; i < (PP); i++)
		{
			V_Col[i] = 0;
			V_Col_Plus[i] = 0;
			V_Col_Minus[i] = 0;
		}
	
		/****************************************************************************************/
		/*概要:ボルテックスの熱散乱・力の初期化													*/
		/*																						*/
		/****************************************************************************************/

		//熱散乱初期化
		gauss->set(0.0, Kt * f0);
		//力の初期化
		calc_f(cd, fc, fc_vvi_p,fc_vvi_m,fc_pin,fc_Lor,fc_tho,vl, ps, side, side_x, side_y);
		//fc0に同じ値を代入
		for (i = 0; i < N * 3; i++)
		{
			fc0[i] = fc[i];
			//printf("%f %f\n",fc0[i],fc[i]);
		}

		/****************************************************************************************/
		/*概要:ピニングサイト配置の初期化														*/
		/*																						*/
		/****************************************************************************************/
#ifdef SET_PINNING_SITE_3
		set_pinning_site_3(ps, side, Rl, Rm, Rs);		
#endif
#ifdef SET_PINNING_SITE_3_2PAIR
		set_pinning_site_3_2Pair(ps, side, Rl, Rm, Rs);
#endif
#ifdef SET_PINNING_SITE_2
		set_pinning_site_2(ps, side, Rl, Rs);
#endif

		

/************************************************************************************************/
/*		概要:座標・速度・力の更新																*/
/*																								*/
/************************************************************************************************/

/********************************/
/*  アニール時間について		*/
/********************************/
		int stepANNEAL = 0;
#ifdef ANNEAL
		stepANNEAL = ANNEAL_TIME;
#endif
		for (step = 1; step <= (total_step+stepANNEAL); step++)
		{
#ifdef LORENTZ
#ifdef CHANGE_KLX
			if ((step % change_pm) == 1) variable = -variable;
			Klx = variable;
			if (step < stepANNEAL) Klx = 0;
#endif


#ifdef CHANGE_KT
			if ((step % change_pm) == 1) Klx = -Klx;
#endif
#endif
#ifdef OUTPUT_CD_VIEW_FILE
			if ((step % 100) == 1) output_file(cd,ps,fc, fc_vvi_p,fc_vvi_m,fc_pin,fc_Lor,fc_tho,side, side_x, side_y,save_num);
#endif
			int vRowCount = -1;
			for (j = 0; j < N * 3; j += 3)
			{
				vl[j + 1] = fc[j + 1] / Eta;
				cd[j + 1] += vl[j + 1] * dt + dt / (2 * Eta) * (fc[j + 1] - fc0[j + 1]);//X軸方向
				vl[j + 2] = fc[j + 2] / Eta;
				cd[j + 2] += vl[j + 2] * dt + dt / (2 * Eta) * (fc[j + 2] - fc0[j + 2]);//Y軸方向

				V += vl[j + 1];
				V_Col[(j / 3) % PP] += vl[j + 1];

				if ((j / 3) % PP == 0) {
					vRowCount++;
				}

				V_Row[vRowCount] += vl[j + 1];

				if (Klx > 0)
				{
					V_plus += vl[j + 1];
				}
				else if(Klx<0) {
					V_minus += vl[j + 1];
				}

				/*
				if (vl[j + 1] > 0)
				{
					V_plus += vl[j + 1];
				}
				else {
					V_minus += vl[j + 1];
				}
				*/
				//printf("%f %f \n", V_plus, V_minus);


				//周期境界条件による折り返し
				if (cd[j + 1] >= side_x) { cd[j + 1] -= side_x; }
				if (cd[j + 1] < 0.0) { cd[j + 1] += side_x; }
				//if (cd[j + 1] < (Dp_S2L / 2.0) - (Rl + 0.1)) { cd[j + 1] += side_x; }
				//if (cd[j + 1] < (Dp_S2L / 2.0) - (Rl+0.1)) { cd[j + 1] += side_x+ (Dp_S2L / 2.0) - (Rl + 0.1); }  //PS_Largeの左端が0未満にはみ出してしまうことを考慮
				if (cd[j + 2] >= side_y) { cd[j + 2] -= side_y; }
				if (cd[j + 2] < 0.0) { cd[j + 2] += side_y; }
			}
			//fc0更新
			for (i = 0; i < N * 3; i++) fc0[i] = fc[i];
			//力の計算
			calc_f(cd, fc, fc_vvi_p,fc_vvi_m,fc_pin,fc_Lor,fc_tho,vl, ps, side, side_x, side_y);
		}

#ifdef VARIABLE
		end = clock();
		double endtime = (max_variable - variable) * (end - start) / (60000 * d_variable);
		printf("計算終了まであと %.3e 分です。\n", endtime);
		printf("%.3e", (end - start) / (60000 * d_variable));
#endif

		double mean_V = 0.0, mean_V_plus = 0.0, mean_V_minus = 0.0;
		mean_V = V / (N * total_step);
		mean_V_plus = V_plus / (N * total_step);
		mean_V_minus = V_minus / (N * total_step);
		//printf("%f %f \n", mean_V_plus, mean_V_minus);

		double real_mean_V = scale_velocity * mean_V;
#ifdef CHANGE_KLX
		FILE* file_1;
		char fn_1[128];
		//errno_t err;
		sprintf_s(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((err = fopen_s(&file_1, fn_1, "at")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		FILE* file_2;
		char fn_2[128];
		//errno_t err;
		sprintf_s(fn_2, "./%s_%d/%s_file_%d/V_Row_Col_%d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, save_num);
		if ((err = fopen_s(&file_2, fn_2, "at")) == NULL)
		{
			fprintf(stderr, "file open error VRowCol %s\n", fn_2);
		}
		
		//fprintf_s(file_1, "%e %e %e %e %e %e %e\n", R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio);
		fprintf_s(file_1, "%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", fabs(Klx) / Kp, mean_V, mean_V_plus, mean_V_minus, mean_V_plus + mean_V_minus,Rl,Rm,Rs,Dp_L2M,Dp_M2S,Dp_S2L);
		fprintf_s(file_2, "%e, %e, %e, %e, %e, %e, %e, %e, %e\n", fabs(Klx) / Kp, mean_V, V_Row[0]/(N*total_step),V_Row[1] / (N * total_step),V_Row[2] / (N * total_step),V_Col[0] / (N * total_step),V_Col[1] / (N * total_step),V_Col[2] / (N * total_step),V_Col[3] / (N * total_step));
		fclose(file_1);
		fclose(file_2);
#endif
#ifdef CHANGE_DP
		FILE* file_1;
		char fn_1[128];
		sprintf(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((file_1 = fopen(fn_1, "a")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		fprintf(file_1, "%e %e %e %e\n", Dp, mean_V, real_Dp, real_mean_V);
		fclose(file_1);
#endif
#ifdef CHANGE_KT
		FILE* file_1;
		char fn_1[128];
		sprintf(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((file_1 = fopen(fn_1, "a")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		fprintf(file_1, "%e %e %e %e\n", Kt, mean_V, temp, real_mean_V);
		fclose(file_1);
#endif
#ifdef CHANGE_TMP
		FILE* file_1;
		char fn_1[128];
		sprintf(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((file_1 = fopen(fn_1, "a")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		fprintf(file_1, "%e %e %e %e\n", period, mean_V, scale_time * period, real_mean_V); //交流の周期を変化させる
		fclose(file_1);
#endif
#ifdef VARIABLE
	}
#endif
}
delete gauss;
return 0;
}


/****************************************************************************************/
/*関数名：output_file																	*/
/*																						*/
/*関数概要：cvdファイル出力																*/
/*																						*/
/*引数		：cd		:現時点での全ボルテックス座標									*/								
/*			  ps		:ピニングサイト半径、座標										*/
/*			  side_x	:ボックス長さX成分												*/
/*			  side_y	:ボックス長さY成分												*/
/*			  save_num	:ファイル記録用番号												*/
/*戻り値	：なし																		*/
/****************************************************************************************/
//2021/10/25	cvdにsidex,sideyを出力させるようにした。
void output_file(double* cd, double*ps ,double* fc,double*fc_vvi_p,double*fc_vvi_m,double*fc_pin,double*fc_Lor,double*fc_tho,double side, double side_x, double side_y,int save_num)
{
	int i;
	FILE* fp;
	//FILE* fp2;
	errno_t err;
	char file_name[128];


	//CD FCプロット用
	char fn_1[128];
	FILE* fp2;
	errno_t err2;

	static int header_flag=0;


	//char file_name2[128];
	static int c = 0;
#ifdef VARIABLE
	sprintf_s(file_name, "./%s_%d/%s_file_%d/variable_%f/cd.%06d.cvd", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, fabs(variable), c++);
	//sprintf_s(file_name2, "./%s_%d/%s_file_%d/variable_%f_ANIME/ANIME_cd%06d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, fabs(variable), c++);
#endif
#ifndef VARIABLE
	sprintf(file_name, "./%s_%d/%s_file_%d/cd.%06d.cvd", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, c++);
	sprintf_s(file_name, "./%s_%d/%s_file_%d/ANIME_cd%06d.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, fabs(variable), c++);
#endif
	sprintf_s(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d_CD_FC.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
	if ((err = fopen_s(&fp, file_name, "wt")) == NULL)
	{
		fprintf(stderr, "file open error %s\n", file_name);
	}
	if ((err2 = fopen_s(&fp2, fn_1, "at")) == NULL)
	{
		fprintf(stderr, "file open error %s\n", fn_1);
	}
	/*if ((err2 = fopen_s(&fp2, file_name2, "wt")) == NULL)
	{
		fprintf(stderr, "file open error %s\n", file_name2);
	}*/
	fprintf(fp, "'box_sx=0 box_sy=0 box_sz=0 box_ex=%f box_ey=%f box_ez=%f box_wt=0.01 \n",//st0=%f st0_pos=(-2.5,-2.5)
		0.0, side_x, side_y);
	for (i = 0; i < N * 3; i += 3)
	{
		//ボルテックス数がピニングサイト数を越えたときはピニングサイト考慮分は0とおいてしまう。ps[N (>=Np)]はアドレスが用意されていないため
		if (i < Np*3) {
			fprintf(fp, "%d 1 %f %f %f %f %f %f %f %f\n", i / 3, cd[i], cd[i + 1], cd[i + 2], ps[i], ps[i + 1], ps[i + 2], side_x, side_y);
		}
		else {
			fprintf(fp, "%d 1 %f %f %f %f %f %f %f %f\n", i / 3, cd[i], cd[i + 1], cd[i + 2], 0.0, 0.0, 0.0, side_x, side_y);
		}
		//

		//fprintf(fp2, "%d %f %f %f %f %f\n", i / 3, cd[i + 1],cd[i + 2],ps[i],ps[i+1],ps[i+2]);
		
		if (header_flag == 0) {	//ヘッダー
			fprintf_s(fp2, "ボルテックス番号, 座標, 合力, VVI作用, 対象VVI, VVI反作用, ピニング力, ローレンツ力, 熱ゆらぎ, 現在速度\n");
			header_flag = 1;
		}
		fprintf_s(fp2, "%d , %f, %f, %f,%f,%f,%f,%f,%f, %f\n", (i/3)+1, cd[i + 1], fc[i + 1], fc_vvi_p[i+1],fc_vvi_m[i],fc_vvi_m[i+1],fc_pin[i+1],fc_Lor[i+1],fc_tho[i+1],fc[i+1]/Eta);
	}
	fclose(fp);

	
	
	fclose(fp2);
}

/****************************************************************************************/
/*関数名：calc_f																		*/
/*																						*/
/*関数概要：ボルテックスにかかる全合力を計算する										*/
/*																						*/
/*引数		：cd	:現時点での全ボルテックス座標										*/
/*			  fc	:現時点での全ボルテックス合力										*/
/*			  fc_***:ボルテックスにかかる力(解析用)										*/										
/*			  vl	:現時点での全ボルテックス速度										*/
/*			  ps	:ピニングサイト半径、座標											*/	
/*			  side	:ボックス長さ														*/	
/*			  side_x:ボックス長さX成分													*/
/*			  side_y:ボックス長さY成分													*/
/*戻り値	：なし																		*/
/****************************************************************************************/

//21.11.01		解析用の変数追加・・・fc_vvi,fc_pin,fc_Lor,fc_tho
//				相互作用、ピニング力、ローレンツ力、熱揺らぎ、それぞれ独立して扱えるようにした。
void calc_f(double* cd, double* fc, double*fc_vvi_p,double*fc_vvi_m,double*fc_pin,double*fc_Lor,double*fc_tho, double* vl, double* ps, double side, double side_x, double side_y)
{
	int pn = N;
	//double a = Dp_S2L;
	int i, j,k;
	double dd[3] = { 0,0,0 };
	double df[3] = { 0,0,0 };

	double ddPin[Np][3];

	double r, rd, tmp;

	double rpin[Np];

	double rNearPinningSite,rdNearPinningSite;

	double sideh = side * 0.5, sideh_x = side_x / 2, sideh_y = side_y / 2;
	double cut_off2 = 100.0;//sideh*sideh;//100.0

	double ddNearPinningSite[3];

	int pinningInflag;
	int pinningEffectiveflag;

	double HealLen=2.0;


	//ピニングポテンシャルの重なり部分を打ち消す際に考慮。
	double SiteSide2SiteSide_SL = Dp_S2L - Rl - Rs;
	double SiteSide2SiteSide_LS = Dp_L2M + Dp_M2S - Rl - Rs;

	



	//初期化
	for (i = 0; i < pn * 3; i++)//i+=3
	{
		fc[i] = 0.0;
		fc_vvi_p[i] = 0.0;
		fc_vvi_m[i] = 0.0;
		fc_pin[i] = 0.0;
		fc_Lor[i] = 0.0;
		fc_tho[i] = 0.0;
	}
	//i番目のボルテックスについて順に力を計算
	for (i = 0; i < pn * 3; i += 3)
	{
		
		//j番目のボルテックスとの相互作用力
		for (j = i + 3; j < pn * 3; j += 3)
		{
			
			dd[1] = cd[i + 1] - cd[j + 1];//x方向　対象のVortexの相互距離
			dd[2] = cd[i + 2] - cd[j + 2];//y方向

			if (dd[1] < -sideh_x) dd[1] += side_x;  //PSL内VortexとPSS内vortexとの距離を最近傍のものにする。(中円をはさまない)
			if (dd[1] > sideh_x) dd[1] -= side_x;  //
			if (dd[2] < -sideh_y) dd[2] += side_y;
			if (dd[2] > sideh_y) dd[2] -= side_y;
			rd = dd[1] * dd[1] + dd[2] * dd[2]; //斜め方向を考慮した相互距離
			r = sqrt(rd);
			if (rd < cut_off2)//fvvカットオフ判定
			{
				tmp = fvv(sqrt(rd))/ r;
				df[1] = tmp * dd[1];
				df[2] = tmp * dd[2];

				fc[i + 1] += df[1];//X作用
				fc_vvi_p[i + 1] += df[1];
				fc[i + 2] += df[2];//Y作用
				fc_vvi_p[i + 2] += df[1];

				fc[j + 1] -= df[1];//X反作用
				fc_vvi_m[i] = (j/3)+1.0;		//どのボルテックスと相互作用しあっているか？(解析用)
				fc_vvi_m[j + 1] -= df[1];
				fc[j + 2] -= df[2];//Y反作用
				fc_vvi_m[j + 1] -= df[2];
			}
		}
#ifdef PINNING_SITE

		//ピニングサイト作用力
		df[1] = 0.0; df[2] = 0.0;
		pinningInflag = 0;
		pinningEffectiveflag = 0;
		for (k = 0; k < Np * 3; k += 3)
		{				

			dd[1] = cd[i + 1] - ps[k + 1];
			dd[2] = cd[i + 2] - ps[k + 2];
			

			if (dd[1] < -sideh_x) dd[1] += side_x;
			if (dd[1] > sideh_x) dd[1] -= side_x;
			if (dd[2] < -sideh_y) dd[2] += side_y;
			if (dd[2] > sideh_y) dd[2] -= side_y;

			rd = dd[1] * dd[1] + dd[2] * dd[2];

			ddPin[k / 3][1] = dd[1];
			ddPin[k / 3][2] = dd[2];



			rpin[k/3] = sqrt(rd);

			if (rpin[k / 3] <= ps[k]) {
				pinningInflag++;
			}
			else if (rpin[k / 3] <= (ps[k] + HealLen)) {
				pinningEffectiveflag++;
			}
		}
		for (k = 0; k < Np * 3; k += 3)
		{
			//カットオフ判定
			/*if (rpin[k/3] < ps[k])
			{
				//pinningInflag++;
				//ボルテックスがピニングサイトの内側にあれば、dfをクリアしてサイト内のみを代入してbreak
				if (r == 0.0)
				{
					tmp = 0.0;
				}
				else {
					tmp = fp_inside(ps[k] - r) / r;//ピニングサイト内側の分
				}
				df[1] = tmp * dd[1];
				df[2] = tmp * dd[2];

				break;
			}*/
			//else 

			//22.01.19 ピニングポテンシャルの重なりを無視できるようなコーディングに着手(未完成)
			
			if((PinningMeshPos_X[k/3][0] < cd[i+1]) && (PinningMeshPos_X[k/3][1] > cd[i+1]) && pinningInflag==0)
			{
				if (rpin[k/3] == 0.0)
				{
					tmp = 0.0;
				}
				else {
					tmp = fp(rpin[k / 3] - ps[k]) / rpin[k / 3];
				}
				
				//サイト外&&カットオフ範囲内ならばdf[]に加算する
				df[1] += tmp * ddPin[k/3][1];
				df[2] += tmp * ddPin[k/3][2];
			}
		}


		fc[i + 1] += df[1];
		fc_pin[i + 1] = df[1];
		fc[i + 2] += df[2];
		fc_pin[i + 2] = df[2];
#endif

		//粘性抵抗っぽいの
		//fc[i+1]-=Eta*vl[i+1];
		//fc[i+2]-=Eta*vl[i+2];
#ifdef LORENTZ
		//ローレンツ力
		fc[i + 1] += Klx * f0; //x方向
		fc_Lor[i+1]= Klx * f0;
		fc[i + 2] += Kly * f0; //y方向
		fc_Lor[i + 2]= Kly * f0;
#endif
#ifdef THERMAL
		//熱振動っぽいの
		fc_tho[i + 1] = gauss->grand();
		fc[i + 1] += fc_tho[i + 1];
		fc_tho[i + 2] = gauss->grand();
		fc[i + 2] += fc_tho[i + 2];
		
#endif
		//奥行方向のクリア
		fc[i] = 0.0;
	}
}


/****************************************************************************************/
/*関数名：set_pinning_site_3															*/
/*																						*/
/*関数概要：指定位置にピニングサイトをセットする										*/
/*																						*/
/*引数		：ps	:ピニングサイト半径、座標											*/
/*			  side	:ボックス長さ														*/
/*			  Rl	:大円の半径															*/
/*			  Rm	:中円の半径															*/
/*			  Rs	:小円の半径															*/
/*戻り値	：なし																		*/
/****************************************************************************************/
void set_pinning_site_3(double* ps, double side, double Rl, double Rm, double Rs)
{
	int i;
	double SiteDis_S2L = Dp_S2L - Rs - Rl;
	double SiteDis_L2M = Dp_L2M - Rl - Rm;
	double SiteDis_M2S = Dp_M2S - Rm - Rs;

	double boxLargeX = (SiteDis_S2L / 2) + (2 * Rl) + (SiteDis_L2M / 2);
	double boxMiddleX = (SiteDis_L2M / 2) + (2 * Rm) + (SiteDis_M2S / 2);
	double boxSmallX = (SiteDis_M2S / 2) + (2 * Rs) + (SiteDis_S2L / 2);
	//初期化
	for (i = 0; i < N * 3; i++)
	{
		ps[i] = 0;
	}
	for (i = 0; i < Np * 3; i += 3)
	{
		//ピニングサイトの大きさps[]をPPごとにスイッチ
		switch ((i / 3) % PP)
		{
		//Rl+0.01としている理由：円の左端=周期境界条件の左端とならないようにするため
		case 0:
			ps[i] = Rl;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_S2L / 2.0;
			ps[i + 1] = Rl+0.01+POS_OFFS;

			//ピニング力が及ぶ範囲を設定
			PinningMeshPos_X[i / 3][0] = 0;
			PinningMeshPos_X[i / 3][1] = boxLargeX;
			
			break;
		case 1:
			ps[i] = Rm;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] =Dp_L2M+ Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M + Rl+0.01+POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxLargeX;
			PinningMeshPos_X[i / 3][1] = boxLargeX + boxMiddleX;

			break;
		case 2:
			ps[i] = Rs;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;

			//ps[i + 1] = Dp_L2M + Dp_M2S + Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M + Dp_M2S + Rl+0.01+POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxLargeX + boxMiddleX;
			PinningMeshPos_X[i / 3][1] = boxLargeX + boxMiddleX + boxSmallX;

			break;
		
		}
		//PSのX座標配置
		//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;//x   //sqrNp*3 = (大中小の数？)×(i,i+1,i+2の数一組)
		//ps[i + 1] = Dp_L * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
		//PSのY座標配置
		ps[i + 2] = Dp * (int)(i / (PP * 3)) + Dp / 2.0;//y

		
		/*FILE* file_2;
		char fn_1[128];
		errno_t err;
		sprintf_s(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d_PS.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((err = fopen_s(&file_2, fn_1, "at")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		//fprintf_s(file_1, "%e %e %e %e %e %e %e\n", R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio);
		fprintf_s(file_2, "%d , %f, %f, %f\n", i,ps[i],ps[i+1],ps[i+2]);
		fclose(file_2);*/
		//cout<<ps[i]<<" "<<ps[i+1]<<" "<<ps[i+2]<<endl;
	}
}


/****************************************************************************************/
/*関数名：set_pinning_site_3_2Pair															*/
/*																						*/
/*関数概要：指定位置にピニングサイトをセットする										*/
/*																						*/
/*引数		：ps	:ピニングサイト半径、座標											*/
/*			  side	:ボックス長さ														*/
/*			  Rl	:大円の半径															*/
/*			  Rm	:中円の半径															*/
/*			  Rs	:小円の半径															*/
/*戻り値	：なし																		*/
/****************************************************************************************/
void set_pinning_site_3_2Pair(double* ps, double side, double Rl, double Rm, double Rs)
{
	int i;
	double SiteDis_S2L = Dp_S2L - Rs - Rl;
	double SiteDis_L2M = Dp_L2M - Rl - Rm;
	double SiteDis_M2S = Dp_M2S - Rm - Rs;

	double boxLargeX = (SiteDis_S2L / 2) + (2 * Rl) + (SiteDis_L2M / 2);
	double boxMiddleX = (SiteDis_L2M / 2) + (2 * Rm) + (SiteDis_M2S / 2);
	double boxSmallX = (SiteDis_M2S / 2) + (2 * Rs) + (SiteDis_S2L / 2);

	//初期化
	for (i = 0; i < N * 3; i++)
	{
		ps[i] = 0;
	}
	for (i = 0; i < Np * 3; i += 3)
	{
		//ピニングサイトの大きさps[]をPPごとにスイッチ
		switch ((i / 3) % PP)
		{
			//Rl+0.01としている理由：円の左端=周期境界条件の左端とならないようにするため
		case 0:
			ps[i] = Rl;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_S2L / 2.0;
			ps[i + 1] = Rl + 0.01 + POS_OFFS;


			//ピニング力が及ぶ範囲を設定
			PinningMeshPos_X[i / 3][0] = 0;
			PinningMeshPos_X[i / 3][1] = boxLargeX;
			break;
		case 1:
			ps[i] = Rm;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] =Dp_L2M+ Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M + Rl + 0.01 + POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxLargeX;
			PinningMeshPos_X[i / 3][1] = boxLargeX + boxMiddleX;

			break;
		case 2:
			ps[i] = Rs;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_L2M + Dp_M2S + Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M + Dp_M2S + Rl + 0.01 + POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxLargeX + boxMiddleX;
			PinningMeshPos_X[i / 3][1] = boxLargeX + boxMiddleX + boxSmallX;

			break;
		case 3:
			ps[i] = Rl;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_S2L / 2.0;
			ps[i + 1] = Dp_S2L+Dp_L2M + Dp_M2S + Rl + 0.01 + POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxLargeX + boxMiddleX + boxSmallX;
			PinningMeshPos_X[i / 3][1] = (2*boxLargeX) + boxMiddleX + boxSmallX;


			break;
		case 4:
			ps[i] = Rm;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] =Dp_L2M+ Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M+ Dp_S2L + Dp_L2M + Dp_M2S + Rl + 0.01 + POS_OFFS;

			PinningMeshPos_X[i / 3][0] = (2 * boxLargeX) + boxMiddleX + boxSmallX;
			PinningMeshPos_X[i / 3][1] = (2 * boxLargeX) + (2 * boxMiddleX) + boxSmallX;

			break;
		case 5:
			ps[i] = Rs;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_L2M + Dp_M2S + Dp_S2L / 2.0;
			ps[i + 1] = Dp_L2M + Dp_M2S + Dp_S2L + Dp_L2M + Dp_M2S + Rl + 0.01 + POS_OFFS;

			PinningMeshPos_X[i / 3][0] = (2 * boxLargeX) + (2 * boxMiddleX) + boxSmallX;
			PinningMeshPos_X[i / 3][1] = (2 * boxLargeX) + (2 * boxMiddleX) + (2 * boxSmallX);

			break;
		}
		//PSのX座標配置
		//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;//x   //sqrNp*3 = (大中小の数？)×(i,i+1,i+2の数一組)
		//ps[i + 1] = Dp_L * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
		//PSのY座標配置
		ps[i + 2] = Dp * (int)(i / (PP * 3)) + Dp / 2.0;//y


		/*FILE* file_2;
		char fn_1[128];
		errno_t err;
		sprintf_s(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d_PS.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((err = fopen_s(&file_2, fn_1, "at")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		//fprintf_s(file_1, "%e %e %e %e %e %e %e\n", R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio);
		fprintf_s(file_2, "%d , %f, %f, %f\n", i,ps[i],ps[i+1],ps[i+2]);
		fclose(file_2);*/
		//cout<<ps[i]<<" "<<ps[i+1]<<" "<<ps[i+2]<<endl;
	}
	//printf("ONESET\n");
}

/****************************************************************************************/
/*関数名：set_pinning_site_2															*/
/*																						*/
/*関数概要：指定位置にピニングサイトをセットする(周期2)									*/
/*																						*/
/*引数		：ps	:ピニングサイト半径、座標											*/
/*			  side	:ボックス長さ														*/
/*			  Rl	:大円の半径															*/
/*			  Rs	:小円の半径															*/
/*戻り値	：なし																		*/
/****************************************************************************************/
void set_pinning_site_2(double* ps, double side, double Rl, double Rs)
{
	int i;
	double SiteDis_S2L = Dp_S2L - Rs - Rl;
	double SiteDis_L2S = Dp_L2M + Dp_M2S - Rl - Rs;

	double boxSmallX = (SiteDis_S2L / 2) + (2 * Rs) + (SiteDis_L2S / 2);
	double boxLargeX = (SiteDis_S2L / 2) + (2 * Rl) + (SiteDis_L2S / 2);

	//初期化
	for (i = 0; i < N * 3; i++)
	{
		ps[i] = 0;
	}

	for (i = 0; i < Np * 3; i += 3)
	{
		//ピニングサイトの大きさps[]をPPごとにスイッチ
		switch ((i / 3) % PP)
		{
			//Rl+0.1としている理由：円の左端=周期境界条件の左端とならないようにするため
		case 0:
			ps[i] = Rs;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] = Dp_S2L / 2.0;
			ps[i + 1] = (SiteDis_L2S/2)+ Rs + 0.01+POS_OFFS;

			//ピニング力が及ぶ範囲を設定
			PinningMeshPos_X[i / 3][0] = 0;
			PinningMeshPos_X[i / 3][1] = boxSmallX;

			break;
		case 1:
			ps[i] = Rl;
			//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
			//ps[i + 1] =Dp_L2M+ Dp_S2L / 2.0;
			ps[i + 1] =(SiteDis_L2S/2)+Rs+ Dp_S2L  + 0.01+POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxSmallX;
			PinningMeshPos_X[i / 3][1] = boxSmallX+boxLargeX;

			break;
		case 2:
			ps[i] = Rs;
				//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;

				//ps[i + 1] = Dp_L2M + Dp_M2S + Dp_S2L / 2.0;
			ps[i + 1] =(SiteDis_L2S/2)+ Dp_S2L+Dp_L2M + Dp_M2S + Rs + 0.01+POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxSmallX + boxLargeX;
			PinningMeshPos_X[i / 3][1] = boxSmallX + boxLargeX + boxSmallX;

			break;
		case 3:
			ps[i] = Rl;
			ps[i + 1] =(SiteDis_L2S/2)+ Dp_S2L+Dp_L2M + Dp_M2S + Rs+Dp_S2L + 0.01+POS_OFFS;

			PinningMeshPos_X[i / 3][0] = boxSmallX + boxLargeX + boxSmallX;
			PinningMeshPos_X[i / 3][1] = boxSmallX + boxLargeX + boxSmallX + boxLargeX;

			break;
				
		}
		//PSのX座標配置
		//ps[i + 1] = Dp * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;//x   //sqrNp*3 = (大中小の数？)×(i,i+1,i+2の数一組)
		//ps[i + 1] = Dp_L * (double)(i % (sqrNp * 3)) / 3.0 + Dp / 2.0;
		//PSのY座標配置
		ps[i + 2] = Dp * (int)(i / (PP * 3)) + Dp / 2.0;//y


		/*FILE* file_2;
		char fn_1[128];
		errno_t err;
		sprintf_s(fn_1, "./%s_%d/%s_file_%d/%s_%d_%d_PS.txt", PROGRAM_NAME, PROGRAM_NUM, PROGRAM_NAME, save_num, PROGRAM_NAME, PROGRAM_NUM, save_num);
		if ((err = fopen_s(&file_2, fn_1, "at")) == NULL)
		{
			fprintf(stderr, "file open error %s\n", fn_1);
		}
		//fprintf_s(file_1, "%e %e %e %e %e %e %e\n", R_standard_min, R_standard_max, dR_standard, Rl_ratio, Rm_ratio, Rs_ratio);
		fprintf_s(file_2, "%d , %f, %f, %f\n", i,ps[i],ps[i+1],ps[i+2]);
		fclose(file_2);*/
		//cout<<ps[i]<<" "<<ps[i+1]<<" "<<ps[i+2]<<endl;
	}
}


/****************************************************************************************/
/*関数名：fvv																			*/
/*																						*/
/*関数概要：Vortex-Vortex-Interaction(VVI)による力を計算する							*/
/*																						*/
/*引数		：ボルテックス位置															*/
/*戻り値	：VVIによってボルテックスにかかる力											*/																					
/****************************************************************************************/
double fvv(const double r)
{
	double f;
	
	f = f0 * pow(2.71828182846, -(r / lambda));
	//f = f0 * gsl_sf_bessel_I1_e(-(r / lambda),);
	return f;
}


/****************************************************************************************/
/*関数名：fp																			*/
/*																						*/
/*関数概要：ピニング力を計算する														*/
/*																						*/
/*引数		：ボルテックス位置															*/
/*戻り値	：ピニング力																*/	
/****************************************************************************************/
//ポテンシャルをtanh(x)の形にしている。その微分がsech^2
double fp(const double r)
{
	double fp;
	fp = -Kp * f0 / gsl_pow_2(cosh(r / Lp));
	return fp;
}

/****************************************************************************************/
/*関数名：fp_step																		*/
/*																						*/
/*関数概要：ピニング力を計算する														*/
/*																						*/
/*引数		：ボルテックス位置															*/
/*戻り値	：ピニング力																*/
/****************************************************************************************/
double fp_step(const double r,const double basis)
{
	double fp;
	if (r > basis) {
		return 0;
	}
	fp = -Kp * f0 / gsl_pow_2(cosh(r / Lp));
	return fp;
}

/****************************************************************************************/
/*関数名：fp_inside																		*/
/*																						*/
/*関数概要：PS内にあるボルテックスのピニング力を計算する								*/
/*																						*/
/*引数		：ボルテックス位置															*/
/*戻り値	：ピニング力																*/																				
/****************************************************************************************/
double fp_inside(double r)
{
	//TODO:Dpカットオフ変更したため、要修正
#ifdef Fp_INSIDE
	if (r > Lp / 2.0)return 0.0;
	return -Kp * f0 * gsl_pow_4(cos(M_PI * r / Lp));
#endif
	return 0.0;
}

#ifdef SQUARE_CD
void set_cd_square_and_vl(double* cd, double* vl, double side, double side_x, double side_y)//初期値の設定
{
	int i;
	for (i = 0; i < N * 3; i += 3)
	{
		cd[i] = 0.0;
		cd[i + 1] = (side_x / NX) * (double)(i % (NX * 3)) / 3.0 + side_x / NX / 2.0;
		cd[i + 2] = (side_y / NY) * (double)(i / (NY * 3)) + side_y / NY / 2.0;

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;
	}
}
#endif
#ifdef TRIANGLE_CD
void set_cd_triangle_and_vl(double* cd, double* vl, double side, double side_x, double side_y)//後で修正
{
	int i;
	for (i = 0; i < N * 3; i += 3)
	{
		cd[i] = 0.0;
		if ((int)(i / (NX * 3)) % 2 == 0) {
			cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp;
		}
		else {
			cd[i + 1] = 0 + ((i % (NX * 3)) / 3) * Dp;
		}
		cd[i + 2] = (sqrt(3.0) * Dp / 2) + (int)(i / (NY * 3)) * sqrt(3.0) * Dp / 2;

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;
	}
}
#endif
#ifdef PS_EDGE_CD
void set_cd_psedge_and_vl(double* cd, double* vl, double side, double side_x, double side_y)//4個以上には対応させていない
{
	int i;
	for (i = 0; i < N * 3; i += 3)
	{
		cd[i] = 0.0;
		cd[i + 2] = Dp / 2 + (int)(i / (NY * 3)) * Dp;
		switch ((i / 3) % PP) {
		case 0:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rl;
			cd[i + 1] = Dp_S2L / 2- Rl;
			//cd[i + 1] = Dp / 2;
			break;
		case 1:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rm;
			cd[i + 1] = Dp_L2M+(Dp_S2L / 2) - Rm;
			//cd[i + 1] = Dp_L2M + (Dp / 2);
			break;
		case 2:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			cd[i + 1] = Dp_L2M + Dp_M2S + (Dp_S2L / 2) - Rs;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp / 2);
			break;
		case 3:
			cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			break;
		case 4:
			cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			break;
		}

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;
	}
}
#endif
#ifdef RANDOM_CD
void set_cd_random_and_vl(double* cd, double* vl, double side, double side_x, double side_y)//後で修正
{
	int i;
	for (i = 0; i < N * 3; i += 3)
	{
		cd[i] = 0.0;
		cd[i + 1] = ((double)rand() / RAND_MAX) * side_x;
		cd[i + 2] = ((double)rand() / RAND_MAX) * side_y;

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;
	}
}
#endif

#ifdef PS_2ND_MF_CD
void set_cd_ps_2ndMF_and_vl(double* cd, double* vl, double side, double side_x, double side_y,double rl,double rm,double rs)
{
	int i;
	int pinning_flag=-1;
	//int pinning_flag = 0;
	for (i = 0; i < N * 3; i += 3)
	{
		if (i % (Np * 3) == 0) pinning_flag = -pinning_flag;
		cd[i] = 0.0;
		cd[i + 2] = Dp / 2 + (int)(i / (NY * 3)) * Dp;
		//cd[i + 2] = Dp / 2 + (int)((i / 3)%PP) * Dp;
		if (cd[i + 2] > side_y) cd[i + 2] - side_y;

		switch ((i / 3) % PP) {
		case 0:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rl;
			//cd[i + 1] = Dp_S2L / 2 -(pinning_flag* Rl/2);
			cd[i + 1] = rl+0.01 - (pinning_flag * rl / 2) + POS_OFFS;
			
			//cd[i + 1] = Dp / 2;
			break;
		case 1:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rm;
			//cd[i + 1] = Dp_L2M + (Dp_S2L / 2) - (pinning_flag * Rm / 2) ;
			cd[i + 1] = Dp_L2M + rl+0.01- (pinning_flag * rm / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + (Dp / 2);
			break;
		case 2:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp_S2L / 2) - (pinning_flag * Rs / 2);
			cd[i + 1] = Dp_L2M + Dp_M2S + rl+0.01 - (pinning_flag * rs / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp / 2);
			break;
		
		}

#ifdef PS_LAST_CD_RAND
		if (i == N * 3 - 3) {
			cd[i+1]= ((double)rand() / RAND_MAX) * side_x;
			cd[i + 2] = ((double)rand() / RAND_MAX) * side_y;
		}
#endif
#ifdef PS_LAST_CD_DISIGNATE
		if (i == N * 3 - 3) {
			cd[i + 1] = (Dp_S2L / 2)+rl+((Dp_L2M-rl-rm)/2);
			//cd[i + 1] = (Dp_S2L / 2) + rl + ((Dp_L2M - rl - rm) / 2);
			cd[i + 2] = Dp*3/4;
		}
#endif 

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;

		
	}
}
#endif

#ifdef PS_2ND_MF_CD_2PAIR
void set_cd_ps_2ndMF_2Pair_and_vl(double* cd, double* vl, double side, double side_x, double side_y, double rl, double rm, double rs)
{
	int i;
	int pinning_flag = -1;
	//int pinning_flag = 0;
	for (i = 0; i < N * 3; i += 3)
	{
		if (i % (Np * 3) == 0) pinning_flag = -pinning_flag;
		cd[i] = 0.0;
		cd[i + 2] = Dp / 2 + (int)(i / (PP * 3)) * Dp;
		//cd[i + 2] = Dp / 2 + (int)((i / 3)%PP) * Dp;
		if (cd[i + 2] > side_y) cd[i + 2] - side_y;

		switch ((i / 3) % PP) {
		case 0:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rl;
			//cd[i + 1] = Dp_S2L / 2 -(pinning_flag* Rl/2);
			cd[i + 1] = rl + 0.01 - (pinning_flag * rl / 2) + POS_OFFS;

			//cd[i + 1] = Dp / 2;
			break;
		case 1:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rm;
			//cd[i + 1] = Dp_L2M + (Dp_S2L / 2) - (pinning_flag * Rm / 2) ;
			cd[i + 1] = Dp_L2M + rl + 0.01 - (pinning_flag * rm / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + (Dp / 2);
			break;
		case 2:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp_S2L / 2) - (pinning_flag * Rs / 2);
			cd[i + 1] = Dp_L2M + Dp_M2S + rl + 0.01 - (pinning_flag * rs / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp / 2);
			break;

		case 3:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rl;
			//cd[i + 1] = Dp_S2L / 2 -(pinning_flag* Rl/2);
			cd[i + 1] = Dp_S2L + Dp_L2M + Dp_M2S + rl + 0.01 - (pinning_flag * rl / 2) + POS_OFFS;

			//cd[i + 1] = Dp / 2;
			break;
		case 4:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rm;
			//cd[i + 1] = Dp_L2M + (Dp_S2L / 2) - (pinning_flag * Rm / 2) ;
			cd[i + 1] = Dp_S2L + Dp_L2M + (2*Dp_M2S) + rl + 0.01 - (pinning_flag * rm / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + (Dp / 2);
			break;
		case 5:
			//cd[i + 1] = Dp / 2 + ((i % (NX * 3)) / 3) * Dp - Rs;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp_S2L / 2) - (pinning_flag * Rs / 2);
			cd[i + 1] = Dp_S2L + (2*Dp_L2M) + (2 * Dp_M2S) + rl + 0.01 - (pinning_flag * rs / 2) + POS_OFFS;
			//cd[i + 1] = Dp_L2M + Dp_M2S + (Dp / 2);
			break;


		}

#ifdef PS_LAST_CD_RAND
		if (i == N * 3 - 3) {
			cd[i + 1] = ((double)rand() / RAND_MAX) * side_x;
			cd[i + 2] = ((double)rand() / RAND_MAX) * side_y;
		}
#endif
#ifdef PS_LAST_CD_DISIGNATE
		if (i == N * 3 - 3) {
			cd[i + 1] = (Dp_S2L / 2) + rl + ((Dp_L2M - rl - rm) / 2);
			//cd[i + 1] = (Dp_S2L / 2) + rl + ((Dp_L2M - rl - rm) / 2);
			cd[i + 2] = Dp * 3 / 4;
		}
#endif 

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;


	}
}
#endif

#ifdef PS_DUO_CD
void set_cd_ps_duo_and_vl(double* cd, double* vl, double side, double side_x, double side_y, double rl, double rm, double rs)
{
	int i;
	int check;

	double SiteDis_L2S= Dp_L2M + Dp_M2S - Rl - Rs;
	double pinning_flag = 0;
	for (i = 0; i < N * 3; i += 3)
	{
		//pinning_flag = 0.0001*i;
		pinning_flag = ((double)rand() / RAND_MAX)-0.5;

		//pinning_flag = -1;
		cd[i] = 0.0;
		cd[i + 2] = Dp / 2 + (int)(i / (PP * 3)) * Dp;
		//cd[i + 2] = Dp / 2 + (int)((i / 3)%PP) * Dp;
		if (cd[i + 2] > side_y) cd[i + 2] - side_y;

		check = (i / 3) % PP;


		//1st以上に入れる際に考慮。対処療法的な処理なので要修正 22/01/08
		if (i >= (Np*3)) {
			if ((i / 3) % 2 == 1) {
				check = 1;
			}
			if ((i / 3) % 2 == 0) {
				check = 3;
			}
			cd[i + 2] = Dp / 2 + (int)(i / (2 * 3)) * Dp;
		}

		switch (check) {
		case 0:
			cd[i + 1] = (SiteDis_L2S / 2)+rs + 0.01 - (pinning_flag * rs / 2);
			break;
		case 1:
			cd[i + 1] = (SiteDis_L2S / 2) + Dp_S2L + rs + 0.01 - (pinning_flag * rl / 2);
			break;

		//以下必要に応じて
		case 2:
			cd[i + 1] = (SiteDis_L2S / 2) + Dp_L2M + Dp_M2S + Dp_S2L + rs + 0.01 - (pinning_flag * rs / 2);
			break;

		case 3:
			cd[i + 1] = (SiteDis_L2S / 2) + Dp_L2M + Dp_M2S + Dp_S2L+ rs+Dp_S2L + 0.01 - (pinning_flag * rl / 2);
			break;
		
		}

#ifdef PS_LAST_CD_RAND
		if (i == N * 3 - 3) {
			cd[i + 1] = ((double)rand() / RAND_MAX) * side_x;
			cd[i + 2] = ((double)rand() / RAND_MAX) * side_y;
		}
#endif
#ifdef PS_LAST_CD_DISIGNATE
		if (i == N * 3 - 3) {
			cd[i + 1] = (Dp_S2L / 2) + rl + ((Dp_L2M - rl - rm) / 2);
			//cd[i + 1] = (Dp_S2L / 2) + rl + ((Dp_L2M - rl - rm) / 2);
			cd[i + 2] = Dp * 3 / 4;
		}
#endif 

		vl[i] = 0.0;
		vl[i + 1] = 0.0;
		vl[i + 2] = 0.0;


	}
}
#endif

double gsl_pow_2(const double x)
{
	return x * x;
}
double gsl_pow_4(const double x)
{
	double x2 = x * x;
	return x2 * x2;
}