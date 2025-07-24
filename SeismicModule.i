%{
    #include "SeismicModule.h"
%}

%include <std_vector.i>

namespace std {
    %template(VectorDouble) std::vector<double>;
}

%include "SeismicModule.h"