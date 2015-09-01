#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

typedef vector<string> stringList;

//[[Rcpp::export]]
stringList string_split(stringList x, char sep, int num){
    int i, j, numChar, start, length, numStrings = x.size();
    stringList result(numStrings);
 
    for (i = 0 ; i < numStrings; i++ ){
        numChar = x[i].length();
        std::vector< int > sepPos(1); // which character is separater
        j = 0;
        for (j = 0 ; j < numChar; j++){
            if (x[i][j] == sep){
                sepPos.push_back(j); 
            }
        }
        sepPos.push_back(numChar); // add last character number
        
        if (num==1){
            start = sepPos[num-1] ;
            length = sepPos[num]-sepPos[num-1];
        }else{
            start = sepPos[num-1] + 1;
            length = sepPos[num]-sepPos[num-1] -1;
        }
        result[i] = x[i].substr(start,length);
    }
    
    return result;
}


//[[Rcpp::export]]
NumericVector removeNA(NumericVector a){
    int vectorLength = a.size();
    NumericVector result(vectorLength);
    for (int i = 0; i < vectorLength ; i++){
        if (NumericVector::is_na(a[i])){
            result[i] = 0;
        }else{
            result[i] = a[i];
        }
    }
    return result;
}
