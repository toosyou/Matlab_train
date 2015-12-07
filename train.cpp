#include <iostream>
#include <cmath>
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

vector<double> gaussian_elimination(const vector< vector<double> >& coefficient,const vector<double>& argument){
    vector<double> ans(coefficient.size(),0.0);
    vector< vector<double> > associated = coefficient;
    for(int i=0;i<coefficient.size();++i)
        associated[i].push_back(argument[i]);

    //gaussian elimination
    int n = coefficient.size();
    for(int i=0; i<n ; i++) {

        // Search for maximum in this column
        double maxEl = std::abs(associated[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (std::abs(associated[k][i]) > maxEl) {
                maxEl = std::abs(associated[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = associated[maxRow][k];
            associated[maxRow][k] = associated[i][k];
            associated[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -associated[k][i]/associated[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    associated[k][j] = 0;
                } else {
                    associated[k][j] += c * associated[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix associated
    for (int i=n-1; i>=0; i--) {
        ans[i] = associated[i][n]/associated[i][i];
        for (int k=i-1;k>=0; k--) {
            associated[k][n] -= associated[k][i] * ans[i];
        }
    }
    return ans;
}

int main(void){
    vector< vector<double> > train_x;
    vector<double> train_t;
    vector< vector< vector<double> > > weight_m2_tmp;
    vector< vector<double> > weight_m2;
    vector<double> tar_m2;
    vector<double> weight_val;

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
        weight_m2[i][6]  *= 2.0;
        weight_m2[i][7]  *= 2.0;
        weight_m2[i][8]  *= 2.0;
        weight_m2[i][10] *= 2.0;
        weight_m2[i][11] *= 2.0;
        weight_m2[i][13] *= 2.0;
    }
    //init tar_m2
    tar_m2.resize(15,0.0);
    //compute tar_m2
    for(int i=0;i<15;++i){
        tar_m2[i] = vector_sum( weight_m2_tmp[0][i] * train_t );
    }

    //compute weight_val by Gaussian elimination
    weight_val = gaussian_elimination(weight_m2,tar_m2);

    cout << "computing done!" <<endl;

    //output all the things
    //output weight_m2
    fstream out_weight_m2("Weight_m2.txt",fstream::out);
    if(!out_weight_m2.is_open()){
        cerr << "***** ERROR : cannot open Weight_m2.txt to write *****" <<endl;
        exit(-1);
    }
    out_weight_m2.precision(12);
    for(int i=0;i<15;++i){
        for(int j=0;j<15;++j){
            out_weight_m2 << fixed << weight_m2[i][j] << "\t";
        }
        out_weight_m2 <<endl;
    }
    out_weight_m2.close();

    //output tar_m2
    fstream out_tar_m2("Tar_m2.txt",fstream::out);
    if(!out_tar_m2.is_open()){
        cerr << "***** ERROR : cannot open Tar_m2.txt to write *****" <<endl;
        exit(-1);
    }
    out_tar_m2.precision(12);
    for(int i=0;i<tar_m2.size();++i){
        out_tar_m2 << fixed << tar_m2[i] << endl;
    }
    out_tar_m2.close();

    //output weight_value.txt
    fstream out_weight_val("weight_val.txt",fstream::out);
    if(!out_weight_val.is_open()){
        cerr << "***** ERROR : cannot open weight_val.txt to write *****" <<endl;
        exit(-1);
    }
    out_weight_val.precision(16);
    for(int i=0;i<weight_val.size();++i){
        out_weight_val << fixed << weight_val[i] << endl;
    }
    out_weight_val.close();

    cout << "output done!" <<endl;
    return 0;
}
