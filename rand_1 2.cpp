//Box-Muller�@�ŗp�����seed���A�萔�l�΂ł͂Ȃ��APC���ԂŌ��܂鐔�l�Ƃ����B

//rand.cpp
//box-muller�@�ɂ�鐳�K���z�����̐���

#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef PI
#define PI 3.141592653589
#endif

//�N���X��
class rand_gauss
{
public:
    //�R���X�g���N�^�G�����o�֐��̏�����
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

    //�����ɂQ�̗����𐶐�����̂ŁA�Ă΂ꂽ���ɃX�g�b�N������΂����f��
    //�X�g�b�N���Ȃ���΁A�V���ɂQ�������Ĉ�f��

    //�����o�֐�
    double grand(void) {
        double temp, r1, r2;

        /*r1 = ((double)rand())/RAND_MAX;
        r2 = ((double)rand())/RAND_MAX;*/
        r1 = ((double)rand() + 1.0) / (RAND_MAX + 2.0);//�͈͂��J���(0,1)�ւƕύX
        r2 = ((double)rand() + 1.0) / (RAND_MAX + 2.0);

        //box-muller�@
        temp = sigma * sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2) + center;
        next = sigma * sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2) + center;
        return temp;
    }
private:
    double center;
    double sigma;
    double next;
};