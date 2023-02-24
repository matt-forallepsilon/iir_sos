#include "iir_sos.h"
#include<array>
#include<vector>
#include<iostream>
#include<numeric>


template<f_type F, c_type C>
void print_filter(Filter<F,C> filter){
    for (auto taps:filter.taps){
        for (auto b:taps.ff_taps)
            cout << b << " ";
        for (auto a:taps.fb_taps)
            cout << a << " ";
        cout << endl;
    }
}

int main()
{
    int n_sections = 3;
    double cutoff = 10;
    double fs = 100;
    double pb_ripple_db = 1;
    double sb_ripple_db = 40;

    Filter<f_type::butter, c_type::lowpass> butter_lowpass(n_sections,cutoff,fs);
    print_filter(butter_lowpass);
    cout << endl;
            
    Filter<f_type::cheby1, c_type::lowpass> cheby1_lowpass(n_sections,cutoff,fs,pb_ripple_db);
    print_filter(cheby1_lowpass);
    cout << endl;

    Filter<f_type::cheby2, c_type::lowpass> cheby2_lowpass(n_sections,cutoff,fs,sb_ripple_db);
    print_filter(cheby2_lowpass);
    cout << endl;

    Filter<f_type::butter, c_type::highpass> butter_highpass(n_sections,cutoff,fs);
    print_filter(butter_highpass);
    cout << endl;
            
    Filter<f_type::cheby1, c_type::highpass> cheby1_highpass(n_sections,cutoff,fs,pb_ripple_db);
    print_filter(cheby1_highpass);
    cout << endl;

    Filter<f_type::cheby2, c_type::highpass> cheby2_highpass(n_sections,cutoff,fs,sb_ripple_db);
    print_filter(cheby2_highpass);
    cout << endl;

    return 0;
}