#pragma once 

#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include "parameters.h"

using namespace std;

class MnsEntry
{
public:
    double val;
    bool is_min;
    double weight;
    int id;
    int index; // position in vector<MnsEntry>

    MnsEntry(){}
    MnsEntry(double _val, bool _is_min, int _id, int _index=-1, double _weight=1):val(_val),
            is_min(_is_min),id(_id),index(_index),weight(_weight){}
};

class MnsRet
{
public:
    double min;
    double max;
    set<int> ids;
    set<int> indexs;
    double weight;

    MnsRet():weight(0){}
};

void reduceVector(vector<MnsEntry> &v, vector<uchar>& status);

int MnS_General(const std::vector<MnsEntry>& _mns_v, std::vector<MnsRet>& mns_ret, int min_inlier=2, bool use_weight=false);