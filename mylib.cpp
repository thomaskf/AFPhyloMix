#include "mylib.h"

string int2str(int i) {
        if (i<0) {
                return "-" + int2str(-i);
        } else if (i<10) {
                return i2s[i];
        } else {
                return int2str(i/10)+i2s[i%10];
        }
}

void tokenizer(string seq, string separators, vector<string>* result) {
    // split the seq into many parts by "separators"
    // the vector<string> *result cannot be NULL
    result->clear();
    int startpos = (int) seq.find_first_not_of(separators);
    while (startpos != (int) string::npos) {
        int endpos = (int) seq.find_first_of(separators, startpos);
        if (endpos != (int) string::npos) {
            result->push_back(seq.substr(startpos, endpos-startpos));
            startpos = (int) seq.find_first_not_of(separators, endpos);
        } else {
            result->push_back(seq.substr(startpos));
            break;
        }
    }
}

void randomizePos(vector<int>& arr, int n) {
        // reorder the items inside the array randomly and only report the top n items
        int i,r,t;
        for (i=0; i<n; i++) {
                r=rand() % (arr.size()-i);
                // swap between pos i and pos r+i
                t = arr[i];
                arr[i] = arr[r+i];
                arr[r+i] = t;
        }
}

void randomizePos(vector<int>& arr) {
        // reorder the items inside the array randomly
	randomizePos(arr, (int)arr.size()-1);
}

void statistics(long double* arr, int n, long double& min, long double& max, long double& avg, long double& median, long double& stdev) {
	// get the minimum, maximum, average, median, and standard deviation for the set of numbers
	vector<long double> a;
	int i;
	for (i=0; i<n; i++)
		a.push_back(arr[i]);
	statistics(a, min, max, avg, median, stdev);
}

void statistics(vector<long double>& arr, long double& min, long double& max, long double& avg, long double& median, long double& stdev) {
	// get the minimum, maximum, average, median, and standard deviation for the set of numbers
	if (arr.size() == 0)
		return;
	int i,n;
	long double t;
	sort(arr.begin(), arr.end());
	long double sum = arr[0];
	min = arr[0];
	max = arr[0];
	n = (int)arr.size();
	for (i=1; i<n; i++) {
		sum += arr[i];
		if (arr[i] < min)
			min = arr[i];
		else if (arr[i] > max)
			max = arr[i];
	}
	avg = sum / n;
	if (n%2 == 1)
		median = arr[(n-1)/2];
	else
		median = (arr[n/2]+arr[n/2-1])/2.0;
	stdev = 0.0;
	for (i=0; i<n; i++) {
		t = arr[i]-avg;
		stdev += (t*t);
	}
	stdev = sqrtl(stdev/(n-1));
}

int getTopK(vector<int>& arr, int topK) {
	// return the value of topK
	vector<int> t;
	int i;
	for (i=0; i<arr.size(); i++)
		t.push_back(arr[i]);
	sort(t.begin(), t.end());
	return t[topK];
}

long subcombin(int n, int k, int dim, long* ans) {
	if (ans[n*dim + k] == 0)
		ans[n*dim + k] = subcombin(n-1, k-1, dim, ans) + subcombin(n-1, k, dim, ans);
	return ans[n*dim + k];
}

long double logfact(int n) {
	
	if (initialCumLogI && n < MAXDIM)
		return cumlogI[n];
		
	long double i = 0.0;
	for (int k=2; k<=n; k++) {
		i += logl((long double)k);
	}
	return i;
}

long double logfact(int n, int k) {
	
	if (initialCumLogI==1 && n < MAXDIM)
		return cumlogI[n] - cumlogI[k-1];
		
	long double i = 0.0;
	for (int m=k; m<=n; m++) {
		i += logl((long double)m);
	}
	return i;
}

long double logCombin(int n, int k) {
	// = log(n!)-log(k!)-log((n-k)!)
	if (n==k)
		return 0.0;
	return logfact(n,k+1)-logfact(n-k);
}

// (1) if p=0 and x>0, return undefine
// (2) if p=1 and x<n, return undefine
long double logBinomial(int x, int n, long double p, bool& isUndefined) {
    if ((p==0.0 && x>0) || (p>=1.0 && x<n)) {
        isUndefined = true;
        return LDBL_MIN_EXP;
    }
    isUndefined = false;
    if (x == 0) {
        return logCombin(n,x) + (n-x)*logl(1.0-p);
    } else if (x == n) {
        return logCombin(n,x) + x*logl(p);
    } else {
        return logCombin(n,x) + x*logl(p) + (n-x)*logl(1.0-p);
    }
}

// (1) if p=0 and x>0, return undefine
// (2) if p=1 and x<n, return undefine
long double logBinomial_noCoeff(int x, int n, long double p, bool& isUndefined) {
    if ((p==0.0 && x>0) || (p>=1.0 && x<n)) {
        isUndefined = true;
        return LDBL_MIN_EXP;
    }
    isUndefined = false;
    if (x == 0) {
        return (n-x)*logl(1.0-p);
    } else if (x == n) {
        return x*logl(p);
    } else {
        return x*logl(p) + (n-x)*logl(1.0-p);
    }
}

/*
long double binomial(int x, int n, long double p) {
    if (n==0)
        return 1.0;
    if (x==0) {
        // i.e. n-x > 0
        if (p>=1.0) { // i.e. (1-p) <= 0
            return 0.0;
        } else {
            return expl(logCombin(n,x) + (n-x)*logl(1.0-p));
        }
    } else if (x==n) {
        // i.e. n-x=0 and x > 0
        if (p==0.0) {
            return 0.0;
        } else {
            return expl(logCombin(n,x) + x*logl(p));
        }
    } else {
        if (p==0.0 || p>=1.0) {
            return 0.0;
        } else {
            return expl(logCombin(n,x) + x*logl(p) + (n-x)*logl(1.0-p));
        }
    }
}
*/

void ratios(vector<long double>& in_logf, vector<long double>& out_ratiof) {
	// given a set of log values,
	// output the ratios between them
	long double logsum;
	long double logratio;
	int i;
	out_ratiof.clear();
	if (in_logf.size() == 0)
		return;
	logsum = in_logf[0];
	for (i=1; i<in_logf.size(); i++) {
		// compute log(x+y), given log(x) and log(y)
		logsum = log_x_plus_y(logsum, in_logf[i]);
	}
	for (i=0; i<in_logf.size(); i++) {
		logratio = in_logf[i] - logsum;
		out_ratiof.push_back(expl(logratio));
	}
}

long combin(int n, int k) {
	// n >= k
	int i;
	if (k>n || n==0)
		return 0;
	if (k==1)
		return n;
	if (k==0 || k==n)
		return 1;
	long *ans = new long[n*k];
	memset(ans, 0, n*k*sizeof(long));
	for (i=0; i<=k; i++)
		ans[i*k+i] = 1;
	for (i=0; i<n; i++)
		ans[i*k] = i+1;
	ans[0] = 1;
	
	long out = subcombin(n-1, k-1, k, ans);
	delete[] ans;
	
	return out;
}

long fact(int n) {
    if (n<=1)
        return 1;
    else
        return n*fact(n-1);
}

void removeSpaces(string& s) {
	// remove all spaces inside the string
	int i,k;
	k=0;
	for (i=0; i<s.length(); i++) {
		if (s[i]!=' ' && s[i]!='-') {
			if (k<i)
				s[k]=s[i];
			k++;
		}
	}
	if (k==0)
		s = "";
	else if (k<s.length())
		s.resize(k);
}

// check whether the string contains space
bool containSpace(string& s) {
	int i;
	for (i=0; i<s.length(); i++)
		if (s[i]==' ' || s[i]=='-')
			return true;
	return false;
}

// is included each other
// 0: no; 1: s1 = s2; 2: s1 is included by s2; 3: s2 is included by s1
int isIncluded(string& s1, string& s2) {
	if (s1 == s2)
		return 1;
	if (s1.length() < s2.length()) {
		if (s2.find(s1) == string::npos)
			return 0;
		else
			return 2;
	} else { // s1.length() > s2.length()
		if (s1.find(s2) == string::npos)
			return 0;
		else
			return 3;
	}
}

// compute the BIC value
long double BIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl ((long double) dataSize);
}

// compute the AIC value
long double AIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * 2.0;
}

// compute the Adjusted BIC value
long double ABIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl (((long double) dataSize + 2.0)/24.0);
}


// integer to string
// minLen : minimum length of the integer.
// If the number of digits < minimum length, then the integer will be displayed with leading zeros, unless the integer is zero
string intToStr(int i, int minLen) {
	char* d2c = (char*) "0123456789"; // the array for digit to char
	if (i<0) {
		return "-" + intToStr(-i, minLen);
	} else if (i<10) {
		if (minLen > 0)
			return string(minLen-1, '0') + string(1,d2c[i]);
		else
			return string(1,d2c[i]);
	} else {
		return intToStr(i/10, minLen-1)+string(1,d2c[i%10]);
	}
}

// double to string
string doublToStr(double d, double decimalPlace) {
    // cout << "[doublToStr enter] d=" << d << " decimalPlace=" << decimalPlace << endl << flush;
    string s;
	if (isnan(d)) {
		s = "N/A";
	} else if (d < 0) {
		s = "-" + doublToStr(-d, decimalPlace);
	} else if (decimalPlace > 0) {
		d = roundNumber(d, decimalPlace);
		s = int2str((int)d) + "." + intToStr((int)round((d - (int)d)*pow(10.0,decimalPlace)), decimalPlace);
	} else {
		s = int2str((int)d);
	}
    // cout << "[doublToStr exit] d=" << d << " decimalPlace=" << decimalPlace << endl << flush;
    return s;
}

// round the double to a certain number of decimal place
double roundNumber(double d, int decimalPlace) {
	return round( d * pow(10.0, decimalPlace) ) / pow(10.0, decimalPlace);
}


// compute the CAIC value
long double CAIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl ((long double) dataSize + 1.0);
}

// compute the AICc value
long double AICc(long double logL, int degFree, int dataSize) {
	return AIC(logL, degFree, dataSize) + (long double) 2.0 * degFree * (degFree + 1) / (dataSize - degFree - 1);
}

// initialize the cumlogI array
void initializeCumLogI() {
	int i;
	cumlogI[1] = 0;
	for (i=2; i<MAXDIM; i++) {
		cumlogI[i] = logl((long double)i) + cumlogI[i-1];
	}
	initialCumLogI = true;
}

// compute log(x+y), given log(x) and log(y)
long double log_x_plus_y(long double logx, long double logy) {
	long double loga, logb;
    
    if (isinf(logx)) {
        return logy;
    }
    if (isinf(logy)) {
        return logx;
    }
	if (logx <= logy) {
		loga = logx;
		logb = logy;
	} else {
		loga = logy;
		logb = logx;
	}
	return logl(expl(loga - logb) + 1) + logb;
}

// compute the log value of dirichlet multinomial distribution
// i.e. log ( prob(beta | alpha) )
// the sizes of alpha and beta are the same
// all elements in alpha has to be greater than one
// need to have invoked initializeCumLogI()
long double logDiriMultDist(vector<int>& alpha, vector<int>& beta) {
	long double result;
	int sum_alpha, sum_beta;
	int i;
	result = 0.0;
	sum_alpha = sum_beta = 0;
	for (i=0; i<alpha.size(); i++) {
		sum_alpha += alpha[i];
		sum_beta += beta[i];
		result += logfact(alpha[i] + beta[i] - 1, alpha[i]);
		result -= logfact(beta[i]);
	}
	result += logfact(sum_beta);
	result -= logfact(sum_alpha + sum_beta - 1, sum_alpha);
	return result;
}

// compute the log value of multinomial distribution
// i.e. log ( prob(freq | p) )
// the sizes of freq and p are the same
long double logMultiDist(vector<int>& freq, vector<double>& p) {
	int sum_freq = 0;
	int i;
	long double result = 0.0;
	for (i=0; i<freq.size(); i++) {
		sum_freq += freq[i];
		result -= logfact(freq[i]);
		result += (long double) freq[i] * logl(p[i]);
	}
	result += logfact(sum_freq);
	return result;
}

// compute the log value of multinomial coefficient
long double logMultiCoeff(vector<int>& freq) {
    int sum_freq = 0;
    int i;
    long double result = 0.0;
    for (i=0; i<freq.size(); i++) {
        sum_freq += freq[i];
        result -= logfact(freq[i]);
    }
    result += logfact(sum_freq);
    return result;
}

// compute the log value of multinomial coefficient
long double logMultiCoeff(int* freq, int n) {
    int sum_freq = 0;
    int i;
    long double result = 0.0;
    for (i=0; i<n; i++) {
        sum_freq += freq[i];
        result -= logfact(freq[i]);
    }
    result += logfact(sum_freq);
    return result;
}

int minInt(int x, int y) {
	if (x<=y)
		return x;
	else
		return y;
}

int maxInt(int x, int y) {
	if (x>=y)
		return x;
	else
		return y;
}

// get maximum of three integers
int maxInt(int x, int y, int z) {
	if (x >= y) {
		if (x >= z) {
			return x;
		} else {
			return z;
		}
	} else {
		if (y >= z) {
			return y;
		} else {
			return z;
		}
	}
}

// get minimum item out of the first K itemss inside doublelist
// numZero - number of items are zeros
// minNoZero - the minimum non-zerio item
void minDouble(vector<double>& doubleList, int& numZero, double& minNoZero, int K) {
    numZero = 0;
    minNoZero = 0.0;
    if (doubleList.size() == 0)
        return;
    int i;
    if (doubleList.size() < K)
        K = doubleList.size();
    for (i=0; i<K; i++) {
        if (doubleList[i] == 0.0)
            numZero++;
        else {
            if (minNoZero == 0.0 || minNoZero > doubleList[i]) {
                minNoZero = doubleList[i];
            }
        }
    }
}

/*
// get minimum from a list of doubles
double minDouble(vector<double>& doubleList, int firstK) {
    if (doubleList.size() == 0)
        return 0.0;
    double min = doubleList[0];
    int i;
    if (doubleList.size() < firstK)
        firstK = doubleList.size();
    for (i=1; i<firstK; i++) {
        if (doubleList[i] < min)
            min = doubleList[i];
    }
    return min;
}
*/

/*
// use continue function to extimate the maximum of a set of long double
// larger value of k_factor, closer to the maximum value
long double maxLongFun(vector<long double>& values, int k_factor) {
    long double log_sum1; // log of sum_i of power (values[i], k_factor-1)
    long double log_sum2; // log of sum_i of power (values[i], k_factor)
    long double log_curr;
    int i;
    if (values.size()==0)
        return 0.0;
    else if (values.size()==1)
        return values[0];
    else {
        log_curr = logl(values[0]);
        log_sum1 = (k_factor-1.0)*log_curr;
        log_sum2 = k_factor*log_curr;
        for (i=1; i<values.size(); i++) {
            log_curr = logl(values[i]);
            log_sum1 = log_x_plus_y(log_sum1, (k_factor-1.0)*log_curr);
            log_sum2 = log_x_plus_y(log_sum2, k_factor*log_curr);
        }
        return expl(log_sum2 - log_sum1);
    }
}
*/

// same as the previous function
// use continue function to extimate the maximum of a set of long double
// larger value of k_factor, closer to the maximum value
// except that all the input values are log values
long double logmaxLongFun(vector<long double>& logvalues, int k_factor) {
    long double log_sum1; // log of sum_i of power (values[i], k_factor-1)
    long double log_sum2; // log of sum_i of power (values[i], k_factor)
    int i;
    if (logvalues.size()==0)
        return 0.0;
    else if (logvalues.size()==1)
        return logvalues[0];
    else {
        log_sum1 = k_factor * logvalues[0];
        log_sum2 = (k_factor+1.0) * logvalues[0];
        for (i=1; i<logvalues.size(); i++) {
            log_sum1 = log_x_plus_y(log_sum1, k_factor * logvalues[i]);
            log_sum2 = log_x_plus_y(log_sum2, (k_factor+1.0) * logvalues[i]);
        }
        return log_sum2 - log_sum1;
    }
}

// get indices with maximum values
void getMaxIndices(vector<long double>& values, vector<int>& maxInd) {
    int i;
    long double max;
    // maxInd.clear();
    if (values.size() == 0)
        return;
    max=values[0];
    // find the max
    for (i=1; i<values.size(); i++) {
        if (values[i] > max)
            max = values[i];
    }
    // find the set of maxIndices
    maxInd.clear();
    for (i=0; i<values.size(); i++) {
        if (fabs(values[i] - max) <= 0.001) {
            maxInd.push_back(i);
        }
    }
    return;
}

// build cigar string
// input: alignS - query alignment
//        alignT - genome alignment
// if the whole read cannot be aligned, return false
bool buildCigar(string& alignS, string& alignT, string& cigarStr) {
	int i, k;
	k=0;
	char pre_state = ' ';
	char cur_state;
	cigarStr = "";
	for (i=0; i<alignS.length(); i++) {
		if (alignS[i] != '-') {
			if (alignT[i] != '-') {
				// match / mismatch
				cur_state = 'M';
			} else if (i==0 || pre_state=='S'){
				// softclipped
				cur_state = 'S';
			} else {
				// insertion
				cur_state = 'I';
			}
		} else {
			// deletion
			cur_state = 'D';
		}
		if (cur_state == pre_state || pre_state == ' ') {
			k++;
		} else {
			cigarStr.append(int2str(k));
			cigarStr.append(1,pre_state);
			k=1;
		}
		pre_state = cur_state;
	}
	if (k>0) {
		if (pre_state == 'S')
			return false;
		if (pre_state == 'I')
			pre_state = 'S'; // softclipped
		cigarStr.append(int2str(k));
		cigarStr.append(1,pre_state);
	}
	return true;
}

// if x1 > x2, swap between them
// else no change
void reorder(int& x1, int& x2) {
	int t;
	if (x1 > x2) {
		t = x1;
		x1 = x2;
		x2 = t;
	}
}

// edit distance between two string
int editDist(string& s1, string& s2) {
	int* s = new int[(s1.length()+1) * (s2.length()+1)];
	// s[i*n + j] = edit distance between s1[1...i] and s2[1...j]
	int i, j, n, v, u, ans;
	n = (int)s2.length() + 1;
	for (i=0; i<=(int)s2.length(); i++)
		s[i] = i;
	for (i=0; i<=(int)s1.length(); i++)
		s[i*n] = i;
	for (i=1; i<=(int)s1.length(); i++) {
		for (j=1; j<=(int)s2.length(); j++) {
			if (s1[i-1] == s2[j-1])
				v = s[(i-1)*n + (j-1)];
			else
				v = s[(i-1)*n + (j-1)]+1;
			u = s[(i-1)*n + j]+1;
			if (u < v)
				v = u;
			u = s[i*n + (j-1)]+1;
			if (u < v)
				v = u;
			s[i*n + j] = v;
		}
	}
	ans = s[s1.length()*n + s2.length()];
	delete[] s;
	return ans;
}

// hamming distance between two strings
// assuming s1.length == s2.length()
int hamDist(string& s1, string& s2) {
	int a = 0;
	int i;
	for (i=0; i<s1.length(); i++) {
		if (s1[i] != s2[i])
			a++;
	}
	return a;
}

#define TOTWEIGHT 10000
// select samples according to weights
void randSelect(vector<double>& weights, vector<int>& selects, int num) {
	double sum = 0;
	double ratio;
	int i,j,k;
	bool isAdded;
	vector<double> cum_weights;
	for (i=0; i<weights.size(); i++) {
		sum += weights[i];
		cum_weights.push_back(sum);
	}
	ratio = (double) TOTWEIGHT / sum;
	for (i=0; i<cum_weights.size(); i++) {
		cum_weights[i] = cum_weights[i] * ratio;
	}
	selects.clear();
	for (i=0; i<num; i++) {
		k = rand() % TOTWEIGHT;
		isAdded = false;
		for (j=0; j<cum_weights.size()-1 && (!isAdded); j++) {
			if (k < cum_weights[j]) {
				selects.push_back(j);
				isAdded = true;
			}
		}
		if (!isAdded)
			selects.push_back((int)cum_weights.size()-1);
	}
}

void timeElapsed(clock_t fr_t, clock_t to_t, string msg) {
	int totSecs;
	int days,hrs,mins,secs;
	bool display;
	totSecs = (int) ((to_t - fr_t)*1.0 / CLOCKS_PER_SEC);
	days = totSecs / 86400;
	totSecs = totSecs % 86400;
	hrs = totSecs / 3600;
	totSecs = totSecs % 3600;
	mins = totSecs / 60;
	secs = totSecs % 60;
	cout << msg << ": ";
	display = false;
	if (days > 0) {
		cout << " " << days << " days";
		display = true;
	}
	if (hrs > 0) {
		cout << " " << hrs << " hours";
		display = true;
	}
	if (mins > 0) {
		cout << " " << mins << " minutes";
		display = true;
	}
	if ((secs > 0)|| (!display))
		cout << " " << secs << " seconds";
	cout << endl;
}

// for Kahan summation algorithm
Summation::Summation() {
    islocked = false;
}

void Summation::reset() {
    s = 0.0L;
    c = 0.0L;
}

void Summation::add(long double input) {
    y = input - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
}

long double Summation::sum() {
    return s;
}

// lock it so that no other function can use it at the same time
void Summation::lock() {
    if (islocked) {
        cerr << "[mylib - Summation] Error! More than one functions are using the same Summation object" << endl;
        exit(1);
    } else {
        islocked = true;
    }
}

// release it so that other function can use it now
void Summation::unlock() {
    islocked = false;
}

// get the double after the position p
// return the double, and the starting and the ending position of the double
double getDouble(string& str, int p, int& fr_p, int& to_p) {
    double d;
    int j = p;
    while (j<str.length() && (!isdigit(str[j])))
        j++;
    fr_p = j;
    while (j<str.length() && (isdigit(str[j]) || str[j]=='.' || str[j]=='-' || str[j]=='e'))
        j++;
    to_p = j;
    if (to_p > fr_p) {
        d = atof(str.substr(fr_p,to_p - fr_p).c_str());
        return d;
    } else {
        cerr << "Error! Double does not exist! " << str.substr(p) << endl;
        exit(1);
    }
}

// generate a set of random numbers
// both fr_value and to_value <= 1.0
void genRandom(double fr_value, double to_value, int n, vector<double>& rand_nums) {
    double d = to_value - fr_value;
    double rand_new;
    int factor = 1000000;
    int d_int = d*factor;
    int i;
    rand_nums.clear();
    // cout << "** d_int=" << d_int << endl;
    if (to_value < fr_value)
        return;
    for (i=0; i<n; i++) {
        rand_new = ((double)(rand() % d_int) / (double)factor) + fr_value;
        // cout << "** " << rand_new << endl;
        rand_nums.push_back(rand_new);
    }
}

void genRandom(double fr_value, double to_value, double sum_value, int n, vector<double>& rand_nums) {
    double sum;
    double ratio;
    int i;
    
    genRandom(fr_value, to_value, n, rand_nums);
    if (rand_nums.size() == 0) {
        return;
    }
    sum = 0.0;
    for (i=0; i<n; i++) {
        sum += rand_nums[i];
    }
    ratio = sum_value / sum;
    for (i=0; i<n; i++) {
        rand_nums[i] = rand_nums[i] * ratio;
    }
}

