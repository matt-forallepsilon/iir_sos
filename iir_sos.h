#ifndef IIR_SOS
#define IIR_SOS

#include<array>
#include<vector>
#include<cmath>
#include <complex> 
#include<iostream>
#include<numeric>
#include<stdexcept>
#include<utility>
#include<limits>


using namespace std;

typedef complex<double> dcomplex;

#define M_I dcomplex(0.0,1.0)

enum class f_type{butter, cheby1, cheby2};
enum class c_type{lowpass,highpass};

struct ZPK{
    vector<dcomplex> zeros;
    vector<dcomplex> poles;
    dcomplex gain = 1.0;
};

struct Biquad{
    array<double,3> ff_taps;
    array<double,3> fb_taps;

    Biquad(const array<double,3>& ff_taps, const array<double,3>& fb_taps){
        this->ff_taps = ff_taps;
        this->fb_taps = fb_taps;
    }
};


void map_cutoff(ZPK& analog_ref, const double cut_freq, const c_type cut_type){
    int N = analog_ref.poles.size();
    if (cut_type == c_type::lowpass){
        if (analog_ref.zeros.size() == 0)
            analog_ref.gain *= pow(cut_freq,N);
        else
            transform(
                analog_ref.zeros.begin(),analog_ref.zeros.end(),analog_ref.zeros.begin(),
                [cut_freq](dcomplex z_k){return cut_freq*z_k;}
            );
        transform(
            analog_ref.poles.begin(),analog_ref.poles.end(),analog_ref.poles.begin(),
            [cut_freq](dcomplex s_k){return cut_freq*s_k;}
        );
    }else if (cut_type == c_type::highpass){
        if (analog_ref.zeros.size() == 0){
            for (auto x:analog_ref.poles)
                analog_ref.gain /= -x;
            analog_ref.zeros = vector<dcomplex>(N,0.0);
        }else{
            for (auto k=0; k<N; k++)
                analog_ref.gain *= analog_ref.zeros[k]/analog_ref.poles[k];
            transform(
                analog_ref.zeros.begin(),analog_ref.zeros.end(),analog_ref.zeros.begin(),
                [cut_freq](dcomplex z_k){return cut_freq/z_k;}
            );
        }
        transform(
            analog_ref.poles.begin(),analog_ref.poles.end(),analog_ref.poles.begin(),
            [cut_freq](dcomplex s_k){return cut_freq/s_k;}
        );
    }else
        throw invalid_argument("cut_freq must be one of c_type::lowpass or c_type::highpass.");

}

void bilinear(ZPK& analog_ref){
    int N = analog_ref.poles.size();
    if (analog_ref.zeros.size() == 0){
        for (auto x:analog_ref.poles)
            analog_ref.gain /= 1.0-x;
        analog_ref.zeros = vector<dcomplex>(N,-1.0); 
    }else{
        for (auto k=0; k<N; k++)
            analog_ref.gain *= (1.0-analog_ref.zeros[k])/(1.0-analog_ref.poles[k]);
        transform(
            analog_ref.zeros.begin(), analog_ref.zeros.end(), analog_ref.zeros.begin(),
            [](dcomplex z_k){return (1.0+z_k)/(1.0-z_k);}
        );
    }
    transform(
        analog_ref.poles.begin(), analog_ref.poles.end(), analog_ref.poles.begin(),
        [](dcomplex s_k){return (1.0+s_k)/(1.0-s_k);}
    );
}

vector<pair<dcomplex,dcomplex>> group_conj(vector<dcomplex> vec){
    vector<pair<dcomplex,dcomplex>> conj_vec;
    while (vec.size() > 0){
        vector<double> magnitudes(vec.size());

        transform(
            vec.begin(), vec.end(), magnitudes.begin(),
            [vec](dcomplex s){return abs(imag(vec[0]*s));}
        );
        
        auto conj_idx = min_element(magnitudes.begin(), magnitudes.end()) - magnitudes.begin();

        if (imag(vec[0])>0)
            conj_vec.push_back(make_pair(vec[0], vec[conj_idx]));
        else
            conj_vec.push_back(make_pair(vec[conj_idx], vec[0]));
        
        vec.erase(vec.begin()+conj_idx);
        vec.erase(vec.begin());
    }
    return conj_vec;
}

vector<dcomplex> conv(const vector<dcomplex>& x, const vector<dcomplex>& y){
    int N = x.size();
    int M = y.size();

    vector<dcomplex> z(N+M-1);

    for (int j=0; j<N; j++)
        for (int k=0; k<M; k++)
            z[j+k] += x[j]*y[k];

    return z;
}

vector<Biquad> zpk_to_sos(const ZPK ref_filt){
    // group conj poles
    auto conj_poles = group_conj(ref_filt.poles);

    // sort poles by closest to the unit circle
    sort(conj_poles.begin(), conj_poles.end(),
        [](pair<dcomplex,dcomplex> a, pair<dcomplex,dcomplex> b){
            return abs(a.first - exp(M_I*arg(a.first))) < 
                abs(b.first - exp(M_I*arg(b.first)));
        }
    );

    // group conj zeros
    auto conj_zeros = group_conj(ref_filt.zeros);

    // sort zeros by closest to poles
    vector<pair<dcomplex,dcomplex>> sorted_zeros;
    for (auto p:conj_poles){
        vector<double> distances(conj_zeros.size());
        transform(
            conj_zeros.begin(), conj_zeros.end(), distances.begin(),
            [p](pair<dcomplex,dcomplex> s){return abs(p.first - s.first);}
        );
        auto idx = min_element(distances.begin(), distances.end()) - distances.begin();
        sorted_zeros.push_back(conj_zeros[idx]);
        conj_zeros.erase(conj_zeros.begin()+idx);
    }

    vector<Biquad> sections;

    for (int i=0; i<sorted_zeros.size(); i++){
        vector<dcomplex> numer = {1.0};
        numer = conv(numer, vector<dcomplex>{1.0, -sorted_zeros[i].first});
        numer = conv(numer, vector<dcomplex>{1.0, -sorted_zeros[i].second});

        array<double,3> real_numer;

        transform(
            numer.begin(),numer.end(),real_numer.begin(),
            [](dcomplex x){return real(x);}
        );

        vector<dcomplex> denom = {1.0};
        denom = conv(denom, vector<dcomplex>{1.0, -conj_poles[i].first});
        denom = conv(denom, vector<dcomplex>{1.0, -conj_poles[i].second});

        array<double,3> real_denom;

        transform(
            denom.begin(),denom.end(),real_denom.begin(),
            [](dcomplex x){return real(x);}
        );
        
        sections.insert(sections.begin(), Biquad(real_numer, real_denom));
    }

    transform(
        sections[0].ff_taps.begin(),sections[0].ff_taps.end(),sections[0].ff_taps.begin(),
        [ref_filt](double x){return real(ref_filt.gain)*x;}
    );

    return sections;
}


template<f_type F> class AnalogReference{};

template<> class AnalogReference<f_type::butter> : public ZPK{
public:
    AnalogReference(const int N)
    {
        this->poles = vector<dcomplex>(N);
        for (auto k=0; k<N; k++)
            this->poles[k] = exp(M_I*M_PI*(2.0*k+1.0+N)/(2.0*N));
    };
};

template<> class AnalogReference<f_type::cheby1> : public ZPK{
public:
    AnalogReference(const int N, const double pass_ripple_db)
    {
        double eps = sqrt(pow(10.0, pass_ripple_db/10.0)-1.0);
        this->gain = 1.0;
        this->poles = vector<dcomplex>(N);
        for (auto k=0; k<N; k++){
            this->poles[k] = M_I*cos((1.0/N)*acos(M_I/eps)+k*M_PI/N);
            this->gain *= this->poles[k];
        }
        this->gain /= sqrt(1.0 + pow(eps, 2.0));
    };
};

template<> class AnalogReference<f_type::cheby2> : public ZPK{
public:
    AnalogReference(const int N, const double stop_ripple_db)
    {
        double eps = 1.0/sqrt(pow(10.0, stop_ripple_db/10.0)-1.0);
        this->gain = 1.0;
        this->poles = vector<dcomplex>(N);
        this->zeros = vector<dcomplex>(N);
        for (auto k=0; k<N; k++){
            this->poles[k] = 1.0/(M_I*cos((1.0/N)*acos(M_I/eps)+k*M_PI/N));
            this->zeros[k] = M_I/cos(M_PI/2.0*(2.0*k+1.0)/N);
            this->gain *= this->poles[k]/this->zeros[k];
        }
    };
};

template<f_type F, c_type C> class Filter{
public:
    vector<Biquad> taps;
    double omega_c;
    int N;

    Filter(const int n_sections, const double cutoff, const double fs)
    {
        this->N = 2*n_sections;
        this->omega_c = tan(M_PI*cutoff/fs);

        AnalogReference<F> ref_filt(this->N);

        map_cutoff(ref_filt, this->omega_c, C);

        bilinear(ref_filt);

        this->taps = zpk_to_sos(ref_filt);
    };
    Filter(const int n_sections, const double cutoff, const double fs, const double ripple_db)
    {
        this->N = 2*n_sections;
        this->omega_c = tan(M_PI*cutoff/fs);

        AnalogReference<F> ref_filt(this->N, ripple_db);

        map_cutoff(ref_filt, this->omega_c, C);

        bilinear(ref_filt);

        this->taps = zpk_to_sos(ref_filt);
    };
};

#endif
