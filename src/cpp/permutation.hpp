#include <map>
#include <cmath>

#define M_LOG2E 1.44269504088896340736 //log2(e)

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
    double myent = 0.0, mysum = 0.0, p = 0.0;
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
    return(myent);
}
