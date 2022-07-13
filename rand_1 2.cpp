//Box-Muller法で用いる種seedを、定数値πではなく、PC時間で決まる数値とした。

//rand.cpp
//box-muller法による正規分布乱数の生成

#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef PI
#define PI 3.141592653589
#endif

//クラス名
class rand_gauss
{
public:
    //コンストラクタ；メンバ関数の初期化
    rand_gauss() {
        srand((unsigned)time(NULL));
        set(0.0, 1.0);
    }
    void set(double center, double sigma) {
        this->center = center;
        this->sigma = sigma;
    }
    double get_center() { return center; }
    double get_sigma() { return sigma; }

    //同時に２つの乱数を生成するので、呼ばれた時にストックがあればそれを吐く
    //ストックがなければ、新たに２個生成して一つ吐く

    //メンバ関数
    double grand(void) {
        double temp, r1, r2;

        /*r1 = ((double)rand())/RAND_MAX;
        r2 = ((double)rand())/RAND_MAX;*/
        r1 = ((double)rand() + 1.0) / (RAND_MAX + 2.0);//範囲を開区間(0,1)へと変更
        r2 = ((double)rand() + 1.0) / (RAND_MAX + 2.0);

        //box-muller法
        temp = sigma * sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2) + center;
        next = sigma * sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2) + center;
        return temp;
    }
private:
    double center;
    double sigma;
    double next;
};