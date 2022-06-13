//
//  main.cpp
//  Gauss-Jordan Elimination
//
//  Created by 김제인 on 2022/06/13.
//

#include <iostream>
#include <vector>
#include <stdlib.h>

const double epsilon = 0.0000000001;

double absolute(double num){
    return num *= num < 0 ? -1 : 1;
}

int32_t non_zero(std::vector<double>& r){
    for(int i = 0; i < r.size(); ++i){
        if(absolute(r[i]) >= epsilon){
            return i;
        }else{
            continue;
        }
    }
    return -1;
}

class matrix{
private:
    std::vector<std::vector<double>> m;
public:
    matrix(){ }~matrix(){ }
    
    std::vector<double> at(int i){
        return this->m[i];
    }
    
    double at(int i, int j){
        return this->m[i][j];
    }
    
    void push_back(std::vector<double> row){
        this->m.push_back(row);
        return;
    }
    
    void mult_row(int row, double k){
        for(auto& i : this->m[row]){
            i *= k;
        }
        this->print();
    }
    
    void mult_row_exchange(int n, int j, double k){
        for(size_t i = 0; i < this->m[0].size(); ++i){
            this->m[j][i] += k * this->m[n][i];
        }
        this->print();
    }
    
    void echelon(){
        for(int i = 0; i < this->m.size()-1; ++i){
            int j = 0;
            if((j = non_zero(this->m[i])) == -1){
                continue;
            }else{
                for(int k = 0; k < this->m.size() - 1 - i; ++k){
                    mult_row_exchange(i, 1 + i + k, -1 * this->m[1 + i + k][j] / this->m[i][j]);
                }
            }
        }
    }
    
    void reduce_echelon(){
        echelon();
        
        for(int i = 0; i < this->m.size(); ++i){
            int j = 0;
            if((j = non_zero(this->m[i])) == -1){
                continue;
            }else{
                mult_row(i, 1 / this->m[i][j]);
            }
        }
        
        for(int i = (int)this->m.size() - 1; i >= 0; --i){
            int j = 0;
            if((j = non_zero(this->m[i])) == -1){
                continue;
            }else{
                for(int k = 0; k < i; ++k){
                    mult_row_exchange(i, i - 1 - k, -1 * this->m[i - 1 - k][j]);
                }
            }
        }
    }
    
    void print(){
        for(int i = 0; i < this->m[0].size() * 7; ++i){
            std::cout << '-';
        }
        std::cout << '\n';
        for(const auto& i : this->m){
            std::cout << "| ";
            for(const auto& j : i){
                std::printf("%5.1lf ", j);
            }
            std::cout << "|\n";
        }
        for(int i = 0; i < this->m[0].size() * 7; ++i){
            std::cout << '-';
        }
        std::cout << '\n';
    }
};

int32_t main(const int32_t argc, const char** argv, const char** env) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    
    // insert code here...
    int size;
    matrix input;
    std::cout << "please enter size of matrix : ";
    std::cin >> size;
    std::cout << "-please enter matrix-\n";
    for(int n = 0; n < size; ++n){
        std::vector<double> row(size);
        for(auto& i : row){
            std::cin >> i;
        }
        input.push_back(row);
    }
    input.reduce_echelon();
    
    return EXIT_SUCCESS;
}
