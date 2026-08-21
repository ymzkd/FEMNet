// SeismicModule_common.i - Common Seismic module definitions (language independent)

// Header is included in fem_common.i

// ResponseMax() の double& 出力引数を各言語ネイティブの出力に変換する。
// (未指定だと double* の不透明ポインタになり、事実上呼び出せない)
// ResponseSpectrum() の同名引数は std::vector<double>& なのでこの %apply は適用されない。
%apply double &OUTPUT { double &max_disp, double &max_accel, double &max_vel };
