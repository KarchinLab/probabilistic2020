#include <map>
#include <cmath>
#include <string>

#define M_LOG2E 1.44269504088896340736L //log2(e)

/* Log base 2 */
inline long double log2(const long double x){
    return  log(x) * M_LOG2E;
}

/* Calculates all position information in one function.
 * ALL OF THE BELOW FUNCTIONS ARE DEPRECATED.
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps positions to number of mutations
 */
std::map<std::string, double> position_info(std::map<int, int> pos_ctr){
    int recurrent_sum = 0, val = 0;
    long double myent_2 = 0.0L, myent_e = 0.0L, mysum = 0.0L, p = 0.0L;
    long double frac_of_uniform_ent = 1.0L, num_pos = 0.0L;
    long double delta_ent = 0.0L;
    std::map<std::string, double> out;
    typedef std::map<int, int>::iterator it_type;

    // count total mutations
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        if (val>1){
            recurrent_sum += val;
        }
        mysum += val;
    }

    // calculate entropy 
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        p = val / mysum;
        myent_2 -= p * log2(p);
        myent_e -= p * log(p);
        num_pos += 1;
    }
    if (num_pos > 1) {
        delta_ent = log(num_pos) - myent_e;
    }
    if (mysum > 1) {
        frac_of_uniform_ent = myent_2 / log2(mysum);
    }

    // put output in a map container
    out["recurrent"] = recurrent_sum;
    out["entropy_fraction"] = frac_of_uniform_ent;
    out["delta_entropy"] = delta_ent;
    return(out);
}


/* Counts the number of recurrent missense mutations.
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps positions to number of mutations
 * 
 * Returns
 * -------
 * mysum : int
 *      number of recurrent mutations
 */
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

/* Calculates the fraction of uniform missense position entropy.
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps positions to number of mutations
 *
 * Returns
 * -------
 * frac_of_uniform_ent : double
 *      Fraction of maximum position entropy defined by
 *      the hypothetical uniform distribution of singleton counts.
 */
double frac_position_entropy(std::map<int, int> pos_ctr) {
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


/* Calculates the difference between the observed missense position 
 * entropy and the entropy under a uniform distribution.
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps positions to number of mutations
 *
 * Returns
 * -------
 * delta_ent : double
 *      computes the difference in missense position entropy between 
 *      the observed and uniform. Natural log is used.
 */
double delta_position_entropy(std::map<int, int> pos_ctr) {
    // define variables
    int val = 0;
    long double myent = 0.0L, mysum = 0.0L, p = 0.0L, num_pos = 0.0L;
    long double delta_ent = 0.0L;
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
        myent -= p * log(p);
        num_pos += 1;
    }
    if (num_pos > 1) {
        delta_ent = log(num_pos) - myent;
    }
    return(delta_ent);
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
