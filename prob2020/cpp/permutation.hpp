#include <map>
#include <cmath>
#include <string>

#define M_LOG2E 1.44269504088896340736L //log2(e)

/* Log base 2 */
inline long double log2(const long double x){
    return  log(x) * M_LOG2E;
}

/* Calculates all position-based statistics in one function.
 * Specifically it calculates:
 * 1. number of recurrent missense mutations
 * 2. fraction of uniform missense entropy
 * 3. delta entropy of missense entropy compared to uniform
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps positions to number of mutations
 *
 * Returns
 * -------
 * out : map<string, double>
 *      STL map contianer containing position statistics
 */
std::map<std::string, double> calc_position_statistics(std::map<int, int> pos_ctr,
                                                       float min_frac=0.02,
                                                       int min_recurrent=2,
                                                       int is_observed=1){
    int recurrent_sum = 0, val = 0, min_frac_thresh = 0;
    // int min_recurrent = 2;  // by default need two recurrent mutations
    long double myent_2 = 0.0L, myent_e = 0.0L, mysum = 0.0L, p = 0.0L;
    long double frac_of_uniform_ent = 1.0L, num_pos = 0.0L;
    long double delta_ent = 0.0L;
    std::map<std::string, double> out;
    typedef std::map<int, int>::iterator it_type;

    // count total mutations
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        mysum += val;
    }

    // set definition of recurrent count number based either on the minimum or 
    // on some percentage of the total missense mutations (specified by min_frac)
    min_frac_thresh = (int) (mysum*min_frac + .99);  // recurrent pos definition using min_frac definition
    min_recurrent = ((min_recurrent>min_frac_thresh) ?  min_recurrent:min_frac_thresh);

    // calculate entropy 
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        // get number of mutations per position
        val = iterator->second;

        // add to recurrent count if defined as recurrently mutated position
        //if (val>=min_recurrent+is_observed){
        if (val>=min_recurrent && is_observed==1){
            recurrent_sum += val;
        } else if(val>=2 && is_observed==0){
            recurrent_sum += val; 
        }

        // update entropy metrics
        //if (val<min_recurrent+is_observed){
        if (val<min_recurrent && is_observed==1){
            p = 1 / mysum;
            for (int i=0; i<val; i++){
                myent_2 -= p * log2(p);
                myent_e -= p * log(p);
                num_pos += 1;
            }
        } else {
            p = val / mysum;
            myent_2 -= p * log2(p);
            myent_e -= p * log(p);
            num_pos += 1;
        }
    }

    // normalize the entropy metrics
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


/* Calculates the effect-based statistics in one function.
 *
 * Parameters
 * ----------
 * pos_ctr : map<int, int>
 *      maps "effects" to number of mutations
 *
 * Returns
 * -------
 * out : map<string, double>
 *      STL map contianer containing position statistics
 */
std::map<std::string, double> calc_effect_statistics(std::map<int, int> pos_ctr,
                                                     float min_frac=0.02,
                                                     int min_recurrent=2,
                                                     int is_observed=1){
    int recurrent_sum = 0, pos = 0, val = 0, min_frac_thresh = 0;
    int num_inactivating = 0;
    // int min_recurrent = 2;  // by default need two recurrent mutations
    long double myent_2 = 0.0L, myent_e = 0.0L, mysum = 0.0L, p = 0.0L;
    long double frac_of_uniform_ent = 1.0L, num_pos = 0.0L;
    long double delta_ent = 0.0L;
    std::map<std::string, double> out;
    typedef std::map<int, int>::iterator it_type;

    // count total mutations
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        val = iterator->second;
        mysum += val;
    }

    // set definition of recurrent count number based either on the minimum or 
    // on some percentage of the total missense mutations (specified by min_frac)
    min_frac_thresh = (int) (mysum*min_frac + .99);  // recurrent pos definition using min_frac defiintion
    min_recurrent = ((min_recurrent>min_frac_thresh) ?  min_recurrent:min_frac_thresh);

    // calculate entropy 
    for(it_type iterator = pos_ctr.begin(); iterator != pos_ctr.end(); iterator++) {
        // get number of mutations per position
        pos = iterator->first;
        val = iterator->second;

        // special case for inactivating mutations
        if (pos==-1){
            num_inactivating = val;
            p = val / mysum;
            myent_2 -= p * log2(p);
            num_pos += 1;
            continue;
        }

        // add to recurrent count if defined as recurrently mutated position
        //if (val>=min_recurrent+is_observed){
        if (val>=min_recurrent && is_observed==1){
            recurrent_sum += val;
        } else if(val>=2 && is_observed==0){
            recurrent_sum += val; 
        }

        // update entropy metrics
        if (val<min_recurrent && is_observed==1){
            p = 1 / mysum;
            for (int i=0; i<val; i++){
                myent_2 -= p * log2(p);
                num_pos += 1;
            }
        } else {
            p = val / mysum;
            myent_2 -= p * log2(p);
            num_pos += 1;
        }
    }

    // normalize the entropy metrics
    if (mysum > 1) {
        frac_of_uniform_ent = myent_2 / log2(mysum);
    }

    // put output in a map container
    out["entropy_fraction"] = frac_of_uniform_ent;
    out["recurrent_sum"] = recurrent_sum;
    out["inactivating_sum"] = num_inactivating;
    return(out);
}
