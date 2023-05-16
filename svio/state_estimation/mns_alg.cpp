#include "mns_alg.h"

void reduceVector(vector<MnsEntry> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

int MnS_General(const std::vector<MnsEntry>& _mns_v, std::vector<MnsRet>& mns_ret, int min_inlier, bool use_weight)
{
    std::vector<MnsEntry> mns_v = _mns_v;

    while (mns_v.size() >= min_inlier)
    {
        MnsRet final_interval, tmp_interval;
        sort(mns_v.begin(), mns_v.end(), [&](const MnsEntry& a, const MnsEntry& b)
            {
                return a.val < b.val;
            });
        for (auto it = mns_v.begin(); it != mns_v.end(); it++)
        {
            if (it->is_min == true)
            {
                tmp_interval.min = it->val;
                tmp_interval.ids.insert(it->id);
                tmp_interval.indexs.insert(it->index);
                tmp_interval.weight += it->weight;
            }
            else
            {
                if ((use_weight == false && tmp_interval.ids.size() > final_interval.ids.size())
                        || (use_weight == true && tmp_interval.weight > final_interval.weight))
                {
                    tmp_interval.max = it->val;    
                    final_interval = tmp_interval;
                }
                tmp_interval.ids.erase(it->id);
                tmp_interval.indexs.erase(it->index);
                tmp_interval.weight -= it->weight;
            }
        }

        if ((use_weight == false && final_interval.ids.size() >= min_inlier)
                || (use_weight == true && final_interval.weight > min_inlier))
        {
            mns_ret.push_back(final_interval);

            vector<uchar> status(mns_v.size(), true);
            for (int i = 0; i < mns_v.size(); i++)
            {
                if (final_interval.ids.find(mns_v[i].id) != final_interval.ids.end())
                    status[i] = false;
            }
            reduceVector(mns_v, status);
        }
        else
            break;
    }

    return mns_ret.size();
}