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


/*
std::vector<int> discrete_uniform_kernel(std::map<int, int> pos_ctr,
                                         int bandwidth){
    int mysum = 0, val, pos, prev_count = 0, first_pos, second_pos;
    std::map<int, int> kernel_counts;
    typedef std::map<int, int>::iterator it_type;

    // get sum of counts
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        mysum += val;
    }

    it_type first_it = pos_ctr.begin();  // first iterator
    rit_type second_it = pos_ctr.rbegin();  // second iterator
    first_pos = first_it->first;
    second_pos = second_it->first;
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        current_pos = iterator->first;
        
        // update iterators and counts
        while(second_pos <= current_pos + bandwidth){
            prev_count += second_it->second;
            second_it++;
            second_pos = second_it->first;
        }
        while(first_pos < current_pos - bandwidth){
            prev_count -= first_it->second;
            first_it++;
            first_pos = first_it->first;
        }

        // save counts within bandwidth 
        kernel_counts.push_back(prev_count)
    }

    return kernel_counts
}
*/
