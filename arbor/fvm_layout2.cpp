#include <memory>

#include "fvm_layout2.hpp"

// Blah. I don't like this interface either.
// Redo:
//
//     'policy' x cable_cell -> bag of mlocations
//     bag of mlocations x morph -> cv_geometry
//
// Policy may want to inspect cable cell for points
// where dynamics change discontinuously.

namespace arb {

struct cv_policy_base {
    virtual locset cv_boundary_points(const cable_cell& cell) const = 0;
    virtual std::unique_ptr<cv_policy_base> clone() const = 0;
    virtual ~cv_policy_base() {}
};

struct cv_policy {
    cv_policy() = delete;
    cv_policy(const cv_policy_base& ref) { // implicit
        policy_ptr = ref.clone();
    }

    cv_policy(cv_policy&&) = default;

    cv_policy(const cv_policy& other):
        policy_ptr(other.policy_ptr->clone()) {}

    cv_policy& operator=(const cv_policy&) = default;
    cv_policy& operator=(cv_policy&&) = default;

    locset cv_boundary_points(const cable_cell& cell) const {
        return policy_ptr->cv_boundary_points(cell);
    }

private:
    std::unique_ptr<cv_policy_base> policy_ptr;
};

// Common flags for CV policies; bitwise composable.
enum class cv_policy_flag {
    none = 0,
    interior_forks = 1<<0,
    single_root_cv = 1<<1
};

struct cv_policy_explicit: cv_policy_base {
    explicit cv_policy_explicit(locset locs):
        locs_(std::move(locs)) {}

    locset cv_boundary_points(const cable_cell&) const override {
        return locs_;
    }

private:
    locset locs_;
};

struct cv_policy_max_extent: cv_policy_base {
    explicit cv_policy_max_extent(double max_extent, cv_policy_flag flags = cv_policy_flag::none)
         max_extent_(max_extent), flags_(flags) {}

    locset cv_boundary_points(const cable_cell&) const override {
        // actual work ...
    }

private:
    cv_policy_flag flags_;
};



}


#if 0
namespace arb {

struct cv_policy_base {
    virtual unsigned append_boundary_points(const em_morphology& em, msize_t branch, std::vector<mlocation>& points) = 0;
};

using cv_policy = std::shared_ptr<cv_policy_base>;

// Common flags for CV policies; bitwise composable.
enum class cv_policy_flag {
    none = 0,
    interior_forks = 1<<0
};

struct cv_policy_max_extent: cv_policy_base {
    cv_policy_max_extent() = delete;
    cv_policy_max_extent(double max_extent, enum cv_policy_flag flags = cv_policy_flag::none):
        max_extent_(max_extent), flags_(flags) {}

    unsigned append_boundary_points(const em_morphology& em, msize_t branch, std::vector<mlocation>& points) override {
        double L = em.branch_length(branch);

        unsigned n = 1;
        if (L>0 && max_extent_>0) {
            n = std::ceil(L/max_extent_);
        }

        if (flags_&cv_policy_flag::interior_forks) {
            for (unsigned i = 0; i<n; ++i) {
                points.push_back(mlocation{branch, (0.5+i)/n});
            }
            return n;
        }
        else {
            for (unsigned i = 0; i<=n; ++i) {
                points.push_back(mlocation{branch, (double)i/n});
            }
            return n+1;
        }
    }

private:
    double max_extent_;
    enum cv_policy_flag flags_;
}

struct cv_policy_fixed_count: cv_policy_base {
    cv_policy_max_extent() = delete;
    cv_policy_max_extent(unsigned cvs_per_branch, enum cv_policy_flag flags = cv_policy_flag::none):
        cvs_per_branch_(max_extent), flags_(flags) {}

    unsigned append_boundary_points(const em_morphology& em, msize_t branch, std::vector<mlocation>& points) override {
        unsigned n = cvs_per_branch_;

        if (flags_&cv_policy_flag::interior_forks) {
            for (unsigned i = 0; i<n; ++i) {
                points.push_back(mlocation{branch, (0.5+i)/n});
            }
            return n;
        }
        else {
            for (unsigned i = 0; i<=n; ++i) {
                points.push_back(mlocation{branch, (double)i/n});
            }
            return n+1;
        }
    }

private:
    unsigned cvs_per_branch_;
    enum cv_policy_flag flags_;
}

struct cv_policy_explicit: cv_policy_base {

}

struct cv_policy_: cv_policy_base {
    cv_policy_max_extent() = delete;
    cv_policy_max_extent(double max_extent, enum cv_policy_flag flags = cv_policy_flag::none):
        max_extent_(max_extent), flags_(flags) {}

    unsigned append_boundary_points(const em_morphology& em, msize_t branch, std::vector<mlocation>& points) override {
        double L = em.branch_length(branch);

        unsigned n = 1;
        if (L>0 && max_extent_>0) {
            n = std::ceil(L/max_extent_);
        }

        if (flags_&cv_policy_flag::interior_forks) {
            for (unsigned i = 0; i<n; ++i) {
                points.push_back(mlocation{branch, (0.5+i)/n);
            }
            return n;
        }
        else {
            for (unsigned i = 0; i<=n; ++i) {
                points.push_back(mlocation{branch, (double)i/n);
            }
            return n+1;
        }
    }

        private:
    double max_extent_;
    enum cv_policy_flag flags_;
}


struct cv_policy {
    enum cv_policy_kind {
        max_length,
        fixed_count
    };


    enum cv_policy_kind policy;
    enum cv_policy_flag flags = cv_policy_flag::none;

    double param = 0;

    cv_policy() = delete;
    cv_policy(enum cv_policy_kind policy, double param, enum cv_policy_flag flags=cv_policy_flag::none):
        policy(policy), param(param), flags(flags)
    {}
};

cv_geometry define_cvs(const em_morphology& em, const locset& boundary_points) {

}

cv_geometry define_cvs(const em_morphology& em, cv_policy pol) {
    auto& morph = em.morph();

    // Walk branches of the morphology in preorder; determine boundary points
    // on each branch via policy.

    std::stack<msize_t> branch_stack;
    branch_stack.push(mnpos);

    while (!branch_stack.empty()) {
        msize_t branch = branch_stack.pop();
        for (auto child: morph.branch_children(branch)) {
            branch_stack.push(child);
        }

        // For now at least, spherical roots are always treated as being
        // entirely within one CV.
        if (branch==0 && morph.spherical_root()) continue;

        switch (pol.policy) {
        case 
    }
}

} // namespace arb
#endif
