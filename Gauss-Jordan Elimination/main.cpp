//
//  main.cpp
//  Gauss-Jordan Elimination
//
//  Created by 김제인 on 2022/06/13.
//

#include <iostream>
#include <vector>
#include <cstdlib>

const double epsilon = 0.0000000001;

using namespace std;

class matrix{
private:
    std::vector<std::vector<double>> m;
    
    inline double absolute(double num){
        return num *= num < 0 ? -1 : 1;
    }
    
    inline int32_t non_zero(std::vector<double>& r){
        for(int i = 0; i < r.size(); ++i){
            if(absolute(r[i]) >= epsilon){
                return i;
            }else{
                continue;
            }
        }
        return -1;
    }
    
    inline int32_t zero(std::vector<double>& row){
        for(int i = 0; i < row.size(); ++i){
            if(row[i] != 0 || row[i] != -0){
                return i;
            }else{
                continue;
            }
        }
        return -1;
    }
    
    inline void row_exchangeB(std::vector<std::pair<int, int>>& rows){
        for(int i = (int)rows.size() - 1; i > 0; --i){
            for(int j = 0; j < i; ++j){
                if(rows[j].second < rows[j + 1].second){
                    auto buffer = rows[j];
                    rows[j] = rows[j + 1];
                    rows[j + 1] = buffer;
                    row_exchange(j, j + 1);
                }
            }
        }
        return;
    }
    
    inline void row_exchangeA(){
        bool ret = false;
        std::vector<std::pair<int, int>> rows(this->m.size());
        for(int i = 0; i < rows.size(); ++i){
            rows[i].second = i;
            rows[i].first = zero(m[i]);
            if(rows[i].first != -1){
                ret = true;
            }else{
                continue;
            }
        }
        if(!ret){
            return;
        }else{
            std::sort(rows.begin(), rows.end(), [](std::pair<int, int> a, std::pair<int, int> b){
                if(a.first == b.first){
                    return a > b;
                }else{
                    return a.first > b.first;
                }
            });
            row_exchangeB(rows);
        }
        return;
    }
public:
    matrix(){ }~matrix(){ }
    
    inline void row_exchange(size_t a, size_t b){
        std::vector<double> temp(this->m[a]);
        for(int i = 0; i < this->m[a].size(); ++i){
            this->m[a][i] = this->m[b][i];
        }
        for(int i = 0; i < this->m[b].size(); ++i){
            this->m[b][i] = temp[i];
        }
        return;
    }
    
    inline std::vector<double> at(int i){
        return this->m[i];
    }
    
    inline double at(int i, int j){
        return this->m[i][j];
    }
    
    inline void push_back(std::vector<double> row){
        this->m.push_back(row);
        return;
    }
    
    inline void mult_row(int row, double k){
        for(auto& i : this->m[row]){
            i *= k;
        }
        this->print();
        return;
    }
    
    inline void mult_row_exchange(int n, int j, double k){
        for(size_t i = 0; i < this->m[0].size(); ++i){
            this->m[j][i] += k * this->m[n][i];
        }
        this->print();
        return;
    }
    
    inline void echelon(){
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
        return;
    }
    
    inline void reduce_echelon(){
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
        return;
    }
    
    inline void print(){
//        row_exchangeA();
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
        return;
    }
};

int32_t main(const int32_t argc, const char** argv, const char** env) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    
    // insert code here...
    int rsize, csize;
    matrix input;
    std::cout << "please enter size of matrix of row : ";
    std::cin >> rsize;
    std::cout << "please enter size of matrix of column : ";
    std::cin >> csize;
    std::cout << "-please enter matrix-\n";
    for(int n = 0; n < csize; ++n){
        std::vector<double> row(rsize);
        for(auto& i : row){
            std::cin >> i;
        }
        input.push_back(row);
    }
    input.reduce_echelon();
    
    return EXIT_SUCCESS;
}
