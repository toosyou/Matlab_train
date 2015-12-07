#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#define num 1369
#define num_t 1369

using namespace std;

vector<double> operator *(const vector<double>& a,const vector<double>& b){
    if(a.size() != b.size())
        return vector<double>(0,0.0);
    vector<double> result(a);
    for(int i=0;i<result.size();++i){
        result[i] = a[i] * b[i];
    }
    return result;
}

double vector_sum(const vector<double>& input){
    double result = 0.0;
    for(int i=0;i<input.size();++i){
        result += input[i];
    }
    return result;
}

int main(void){
    vector< vector<double> > train_x;
    vector<double> train_t;
    vector< vector< vector<double> > > weight_m2_tmp;
    vector< vector<double> > weight_m2;
    vector<double> tar_m2;

    //load train_x
    fstream in_train_x("train.txt",fstream::in);
    if(!in_train_x.is_open()){
        cerr << "******* ERROR : cannot open train.txt *******" <<endl;
        exit(-1);
    }
    train_x.resize(4);
    for(int i=0;i<4;++i)
        train_x[i].resize(num,0.0);
    for(int i=0;i<num;++i){
        for(int j=0;j<4;++j){
            in_train_x >>  train_x[j][i];
        }
    }
    in_train_x.close();

    //load train_t
    fstream in_train_t("train_target.txt",fstream::in);
    if(!in_train_t.is_open()){
        cerr << "******* ERROR : cannot open train_target.txt *******" <<endl;
        exit(-1);
    }
    train_t.resize(num_t,0.0);
    for(int i=0;i<num_t;++i){
        in_train_t >> train_t[i];
    }
    in_train_t.close();

    //init weight_m2 & weight_m2_tmp
    weight_m2.resize(15);
    for(int i=0;i<15;++i)
        weight_m2[i].resize(15,0);
    weight_m2_tmp.resize(15);
    for(int i=0;i<15;++i){
        weight_m2_tmp[i].resize(15);
        for(int j=0;j<15;++j){
            weight_m2_tmp[i][j].resize(num,1.0);
        }
    }

    //compute weight_m2_tmp
    //border first
    for(int i=0;i<4;++i){
        weight_m2_tmp[0][1+i] = train_x[i];
        weight_m2_tmp[1+i][0] = train_x[i];
    }
    for(int i=0,shift=0;i<4;++i){
        for(int j=i;j<4;++j){
            vector<double> tmp = train_x[i] * train_x[j];
            weight_m2_tmp[0][5+shift] = tmp;
            weight_m2_tmp[5+shift][0] = tmp;
            shift++;
        }
    }
    //content
    for(int i=1;i<15;++i){
        for(int j=1;j<15;++j){
            weight_m2_tmp[i][j] = weight_m2_tmp[i][0] * weight_m2_tmp[0][j];
        }
    }
    //compute weight_m2. weight_m2[i][j] = sum( weight_m2_tmp[i][j] )
    for(int i=0;i<15;++i){
        for(int j=0;j<15;++j){
            weight_m2[i][j] = vector_sum(weight_m2_tmp[i][j]);
        }
    }
    //do the multiply-2-things
    for(int i=0;i<15;++i){
        weight_m2[i][7] *= 2;
        weight_m2[i][8] *= 2;
        weight_m2[i][9] *= 2;
        weight_m2[i][11] *= 2;
        weight_m2[i][12] *= 2;
        weight_m2[i][14] *= 2;
    }
    //init tar_m2
    tar_m2.resize(15,0.0);
    //compute tar_m2
    for(int i=0;i<15;++i){
        tar_m2[i] = vector_sum( weight_m2_tmp[0][i] * train_t );
    }


    return 0;
}
