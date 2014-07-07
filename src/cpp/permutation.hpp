#include <map>
#include <cmath>

#define M_LOG2E 1.44269504088896340736L //log2(e)

inline long double log2(const long double x){
    return  log(x) * M_LOG2E;
}

int recurrent_sum(std::map<int, int> pos_ctr) {
    int mysum = 0, val = 0;
    typedef std::map<int, int>::iterator it_type;
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        if (val>1){
            mysum += val;
        }
    }
    return(mysum);
}

double position_entropy(std::map<int, int> pos_ctr) {
    // define variables
    int val = 0;
    long double myent = 0.0L, mysum = 0.0L, p = 0.0L;
    long double frac_of_uniform_ent = 1.0L;
    typedef std::map<int, int>::iterator it_type;

    // calculate sum
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        mysum += val;
    }

    // calculate entropy 
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        p = val / mysum;
        myent -= p * log2(p);
    }
    if (mysum > 1) {
        frac_of_uniform_ent = myent / log2(mysum);
    }
    return(frac_of_uniform_ent);
}
