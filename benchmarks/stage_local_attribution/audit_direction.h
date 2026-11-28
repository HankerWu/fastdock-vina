// Independent diagnostic of a single-ligand rigid receptor research panel.
// Checks the actual candidate tangent against translation and angle partitions
// obtained directly from the starting conformation and raw BFGS direction.
#pragma once
#include <fstream>
#include <iomanip>
#include <cstdlib>

struct direction_audit {
    unsigned long unchanged=0, preserve=0, fallback=0, zero_translation_work=0;
    unsigned long nonnegative_before=0, nonnegative_after=0;
    fl max_work_error=0;
    ~direction_audit() {
        std::ofstream file(std::getenv("FD_AUDIT_DIRECTION"));
        file << "small_or_pure,preserve,full_scale_fallback,zero_trans_work_when_capped,"
                "nonnegative_before,nonnegative_after,max_relative_slope_identity_error\n"
             << unchanged << ',' << preserve << ',' << fallback << ',' << zero_translation_work
             << ',' << nonnegative_before << ',' << nonnegative_after << ','
             << std::setprecision(17) << max_work_error << '\n';
    }
    void check(const conf& c,const flv& g,const flv& raw,const flv& proposal) {
        VINA_CHECK(c.ligands.size()==1 && c.flex.empty() && raw.size()==6+c.ligands.front().torsions.size());
        const double rotation=std::hypot(raw[3],std::hypot(raw[4],raw[5]));
        double max_angle=rotation;
        for (sz j=6; j<raw.size(); ++j) max_angle=std::max(std::abs(raw[j]),max_angle);
        const double k=max_angle>pi ? pi/max_angle : 1;
        const double t=raw[0]*g[0]+raw[1]*g[1]+raw[2]*g[2];
        if (k==1) ++unchanged;
        else if (t<=0) {
            ++preserve;
            if (t==0) ++zero_translation_work;
        } else ++fallback;
        double initial=0, final=0, l1=0;
        for (sz j=0; j<g.size(); ++j) {
            const double expected=raw[j]*(j>=3 || t>0 ? k : 1);
            VINA_CHECK(std::isfinite(proposal[j]) &&
                std::abs(proposal[j]-expected)<=8*std::numeric_limits<double>::epsilon()*std::max(fl(1),std::abs(expected)));
            initial+=g[j]*raw[j];
            final+=g[j]*proposal[j];
            l1+=std::abs(g[j]*proposal[j]);
        }
        VINA_CHECK(std::isfinite(initial) && std::isfinite(final));
        nonnegative_before+=initial>=0;
        nonnegative_after+=final>=0;
        const fl expected=k*initial+(t<=0 ? (1-k)*t : 0);
        const double discrepancy=std::abs(expected-final)/std::max(fl(1),l1);
        max_work_error=std::max(max_work_error,discrepancy);
        VINA_CHECK(discrepancy<1e-12 && (initial>=0 || final<0));
    }
};
