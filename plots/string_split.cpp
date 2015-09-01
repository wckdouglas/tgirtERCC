#include <Rcpp.h>
#include <sstream>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace std;

typedef vector<string> stringList;
typedef vector<int> numList;

//split function to split line with desired deliminator
stringList split(const string &s, char delim)
{
    stringList result;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim))
    {
        result.push_back(item);
    }
    return result;
}

//[[Rcpp::export]]
stringList string_split(stringList x, string sep, int start, int frag)
{
        int size = x.size();
        stringList result(size);
        char delim = sep[0];
		for (int i = 0; i < size; i++)
        {
                stringList splittedString = split(x[i],delim);
                int splittedSize = splittedString.size();
				/*
                if (splittedSize < start + frag)
                {
                        Rcpp::stop("Wrong end value\n");
                }
				*/
                string str = splittedString[start-1];
                for (int j = start ; j < start + frag ; j++)
                {
                        str = str + sep + splittedString[j];
                }
                result[i] = str;
        }
    return result;
}

//[[Rcpp::export]]
stringList changeDelim(stringList x, char sep, string delim)
{
        int size = x.size();
        stringList result(size);
        for (int i = 0 ; i < size ; i ++)
        {
                stringList splittedString = split(x[i],sep);
                string str = splittedString[0];
                for (int j = 1 ; j < splittedString.size(); j++)
                {
                        str = str + delim + splittedString[j];
                }
                result[i] = str;
        }
        return result;
}
