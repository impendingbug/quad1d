/**
 *  \file  quad1d.hpp
 *  \brief  Wrapper for cr_quad1d.h and some of GSL quadratures
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  12-12-15
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#ifndef quad1d_hpp_INCLUDED
#define quad1d_hpp_INCLUDED

#include "cr_quad1d.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <boost/current_function.hpp>

#include <cstddef>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <limits>
#include <type_traits>

#define ERROR_STR(str) (std::string(BOOST_CURRENT_FUNCTION) + " : " + str)

namespace quad1d {

template<typename Res, typename Func>
class Bone {
public:
    Bone(): merr_h(gsl_set_error_handler_off()) { }
    virtual ~Bone() = default;

    virtual void
        integrate(const Func& func, double a, double b, Res& res, Res& abserr, std::size_t) = 0;
    virtual Res
        integrate(const Func& func, double a, double b, Res& abserr, std::size_t) = 0;

protected:
    gsl_error_handler_t *merr_h {nullptr};
    const Func *mpfunc {nullptr};
};

//**************************************************************************************

//Available Gauss-Kronrod rules
enum class GK_rule: int { GK15 = 1, GK21 = 2, GK31 = 3, GK41 = 4, GK51 = 5, GK61 = 6 };

template<typename Res, typename Func>
class Bone_ag: public Bone<Res, Func> {
    using base = Bone<Res, Func>;
protected:
    struct Neg_inf { explicit operator double() { return std::numeric_limits<double>::quiet_NaN(); } };
    struct Pos_inf { explicit operator double() { return std::numeric_limits<double>::quiet_NaN(); } };
public:
    Bone_ag() = delete;
    Bone_ag(Res epsrel, Res epsabs, std::size_t max_limit, GK_rule rule):
        base(), mepsrel(epsrel), mepsabs(epsabs), mmax_limit(max_limit), mrule(rule) { }

    template<
        typename Ta, typename Tb,
        typename = typename std::enable_if<
                       (std::is_arithmetic<Ta>::value && std::is_arithmetic<Tb>::value) ||
                       (std::is_same<Neg_inf, Ta>::value && std::is_arithmetic<Tb>::value) ||
                       (std::is_arithmetic<Ta>::value && std::is_same<Pos_inf, Tb>::value) || 
                       (std::is_same<Neg_inf, Ta>::value && std::is_same<Pos_inf, Tb>::value)
                   >::type
    >
    void
    integrate(const Func& func, Ta a, Tb b, Res& res, Res& abserr, std::size_t limit = 2000) {
        this->merr_h = gsl_set_error_handler_off();
        this->mpfunc = &func;
        auto status = integrate_impl(static_cast<double>(a), std::is_arithmetic<Ta>(), static_cast<double>(b),
                                     std::is_arithmetic<Tb>(), res, abserr, limit);
        gsl_set_error_handler(this->merr_h);
        this->handle_err_status(status);  //throw exception (or not) based on status
    }

    template<typename Ta, typename Tb>
    Res
    integrate(const Func& func, Ta a, Tb b, Res& abserr, std::size_t limit = 2000) {
        Res res;
        integrate(func, a, b, res, abserr, limit);
        return res;
    }

    void set_rule(GK_rule rule) noexcept { mrule = rule; };
    void set_acc(Res epsrel, Res epsabs) noexcept { mepsrel = epsrel; mepsabs = epsabs; }
    Neg_inf neg_inf() const noexcept { return Neg_inf{}; }
    Pos_inf pos_inf() const noexcept { return Pos_inf{}; }

protected:
    Res mepsrel;
    Res mepsabs;
    std::size_t mmax_limit;
    GK_rule mrule;

    static Res
    func_wrapper(double x, void *pq) {
        auto pquad = static_cast<Bone_ag *>(pq);
        return (*pquad->mpfunc)(x);
    }

    void handle_err_status(int status) const {
        auto str_err = gsl_strerror(status);
        switch (status) {
            case GSL_SUCCESS: break;
            case GSL_EINVAL: throw std::invalid_argument(ERROR_STR("iteration limit exceeds available workspace"));
            case GSL_EBADTOL: throw std::invalid_argument(ERROR_STR(str_err));
            case GSL_ESING: throw std::domain_error(ERROR_STR(str_err));
            case GSL_EMAXITER: throw std::runtime_error(ERROR_STR(str_err));
            case GSL_EROUND: throw std::runtime_error(ERROR_STR("roundoff error prevents requested tolerance from "
                                                                "being achieved"));
            case GSL_EDIVERGE: throw std::runtime_error(ERROR_STR("integrand is divergent, or slowly convergent"));
            default: throw std::logic_error(ERROR_STR(str_err));  //never thrown
        }
    }

    virtual int
    integrate_impl(double a, std::true_type, double b, std::true_type, Res& res, Res& abserr,
                   std::size_t limit) { return 0; }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::false_type, Res& res, Res& abserr,
                   std::size_t limit) { return 0; }
    virtual int
    integrate_impl(double a, std::true_type, double b, std::false_type, Res& res, Res& abserr,
                   std::size_t limit) { return 0; }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::true_type, Res& res, Res& abserr,
                   std::size_t limit) { return 0; }
};

template<typename Func = std::complex<double>(double)>
class Cag: public Bone_ag<std::complex<double>, Func> {
protected:
    using cx_double = std::complex<double>;
    using base = Bone_ag<cx_double, Func>;
public:
    explicit Cag(cx_double epsrel = {1.E-10, 1.E-10}, cx_double epsabs = 0., std::size_t max_limit = 2000,
                 GK_rule rule = GK_rule::GK31):
        base(epsrel, epsabs, max_limit, rule),
        mw_re(qag_workspace_alloc(max_limit), qag_workspace_free),
        mw_im(qag_workspace_alloc(max_limit), qag_workspace_free)
    {
        if (mw_re.get() == nullptr || mw_im.get() == nullptr) {
            gsl_set_error_handler(this->merr_h);
            throw std::bad_alloc{};
        }
        mf.function = &this->func_wrapper;
        mf.params = this;
        gsl_set_error_handler(this->merr_h);
    }

    Cag(double epsrel, double epsabs, std::size_t max_limit = 2000, GK_rule rule = GK_rule::GK31): 
        Cag( cx_double(epsrel, epsrel), cx_double(epsabs, epsabs), max_limit, rule ) { }

    Cag(const Cag&) = delete;
    Cag(Cag&&) = default;
    Cag& operator=(const Cag&) = delete;
    Cag& operator=(Cag&&) = default;

    using base::integrate;

    virtual void
    integrate(const Func& func, double a, double b, cx_double& res, cx_double& abserr, std::size_t limit = 2000) 
    override {
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
    }

    virtual cx_double
    integrate(const Func& func, double a, double b, cx_double& abserr, std::size_t limit = 2000) override {
        cx_double res;
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
        return res;
    }

    using base::set_acc;
    void set_acc(double epsrel, double epsabs) { base::set_acc(cx_double(epsrel, epsrel), cx_double(epsabs, epsabs)); }

protected:
    using workspace_t = std::unique_ptr<qag_workspace, decltype(&qag_workspace_free)>;
    workspace_t mw_re;
    workspace_t mw_im;
    cr_function mf;

    virtual int
    integrate_impl(double a, std::true_type, double b, std::true_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return cr_gsl_integration_qag(&this->mf, a, b, this->mepsabs, this->mepsrel, limit,
                                      static_cast<int>(this->mrule),
                                      this->mw_re.get(), this->mw_im.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::false_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return cr_gsl_integration_qagi(&this->mf, this->mepsabs, this->mepsrel, limit,
                                       static_cast<int>(this->mrule),
                                       this->mw_re.get(), this->mw_im.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::true_type, double b, std::false_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return cr_gsl_integration_qagiu(&this->mf, a, this->mepsabs, this->mepsrel, limit,
                                        static_cast<int>(this->mrule),
                                        this->mw_re.get(), this->mw_im.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::true_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return cr_gsl_integration_qagil(&this->mf, b, this->mepsabs, this->mepsrel, limit,
                                        static_cast<int>(this->mrule),
                                        this->mw_re.get(), this->mw_im.get(), &res, &abserr);
    }
};

template<typename Func = double(double)>
class Rag_gsl: public Bone_ag<double, Func> {
protected:
    using base = Bone_ag<double, Func>;
public:
    explicit Rag_gsl(double epsrel = 1.E-10, double epsabs = 0., std::size_t max_limit = 2000,
                     GK_rule rule = GK_rule::GK31):
        base(epsrel, epsabs, max_limit, rule),
        mw(gsl_integration_workspace_alloc(max_limit), gsl_integration_workspace_free)
    {
        if (mw.get() == nullptr) {
            gsl_set_error_handler(this->merr_h);
            throw std::bad_alloc{};
        }
        mf.function = &this->func_wrapper;
        mf.params = this;
        gsl_set_error_handler(this->merr_h);
    }

    Rag_gsl(const Rag_gsl&) = delete;
    Rag_gsl(Rag_gsl&&) = default;
    Rag_gsl& operator=(const Rag_gsl&) = delete;
    Rag_gsl& operator=(Rag_gsl&&) = default;

    using base::integrate;

    virtual void
    integrate(const Func& func, double a, double b, double& res, double& abserr, std::size_t limit = 2000) 
    override {
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
    }

    virtual double
    integrate(const Func& func, double a, double b, double& abserr, std::size_t limit = 2000) override {
        double res;
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
        return res;
    }

    void
    integrate_sing(const Func& func, double a, double b, double &res, double& abserr, std::size_t limit = 2000) {
        this->merr_h = gsl_set_error_handler_off();
        this->mpfunc = &func;
        auto status = gsl_integration_qags(&this->mf, a, b, this->mepsabs, this->mepsrel, limit, this->mw.get(), &res,
                                           &abserr);
        gsl_set_error_handler(this->merr_h);
        this->handle_err_status(status);  //throw exception (or not) based on status
    }

    double
    integrate_sing(const Func& func, double a, double b, double& abserr, std::size_t limit = 2000) {
        double res;
        integrate_sing(func, a, b, res, abserr, limit);
        return res;
    }
    
    using base::set_acc;

protected:
    using workspace_t = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;
    workspace_t mw;
    gsl_function mf;

    virtual int
    integrate_impl(double a, std::true_type, double b, std::true_type,
                   double& res, double& abserr, std::size_t limit) override {
        return gsl_integration_qag(&this->mf, a, b, this->mepsabs, this->mepsrel, limit,
                                   static_cast<int>(this->mrule), this->mw.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::false_type,
                   double& res, double& abserr, std::size_t limit) override {
        return gsl_integration_qagi(&this->mf, this->mepsabs, this->mepsrel, limit, this->mw.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::true_type, double b, std::false_type,
                   double& res, double& abserr, std::size_t limit) override {
        return gsl_integration_qagiu(&this->mf, a, this->mepsabs, this->mepsrel, limit, this->mw.get(), &res, &abserr);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::true_type,
                   double& res, double& abserr, std::size_t limit) override {
        return gsl_integration_qagil(&this->mf, b, this->mepsabs, this->mepsrel, limit, this->mw.get(), &res, &abserr);
    }
};

template<typename Func = std::complex<double>(double)>
class Cag_gsl: public Bone_ag<std::complex<double>, Func> {
protected:
    using cx_double = std::complex<double>;
    using base = Bone_ag<cx_double, Func>;
public:
    explicit Cag_gsl(cx_double epsrel = {1.E-10, 1.E-10}, cx_double epsabs = 0., std::size_t max_limit = 2000,
                     GK_rule rule = GK_rule::GK31):
        base(epsrel, epsabs, max_limit, rule), mquad(epsrel.real(), epsabs.real(), max_limit, rule)
    {
        gsl_set_error_handler(this->merr_h);
    }

    Cag_gsl(double epsrel, double epsabs, std::size_t max_limit = 2000, GK_rule rule = GK_rule::GK31): 
        Cag_gsl( cx_double(epsrel, epsrel), cx_double(epsabs, epsabs), max_limit, rule ) { }

    Cag_gsl(const Cag_gsl&) = delete;
    Cag_gsl(Cag_gsl&&) = default;
    Cag_gsl& operator=(const Cag_gsl&) = delete;
    Cag_gsl& operator=(Cag_gsl&&) = default;

    using base::integrate;

    virtual void
    integrate(const Func& func, double a, double b, cx_double& res, cx_double& abserr, std::size_t limit = 2000) 
    override {
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
    }

    virtual cx_double
    integrate(const Func& func, double a, double b, cx_double& abserr, std::size_t limit = 2000) override {
        cx_double res;
        base::template integrate<double, double>(func, a, b, res, abserr, limit);
        return res;
    }

    void
    integrate_sing(const Func& func, double a, double b, cx_double &res, cx_double& abserr, std::size_t limit = 2000) {
        Re_t<Func> re_func(func);
        Im_t<Func> im_func(func);
        double re, im, err_re, err_im;

        mquad.set_acc(this->mepsrel.real(), this->mepsabs.real());
        try {
            mquad.integrate_sing(re_func, a, b, re, err_re, limit);
        } catch (std::runtime_error& e) {
            mquad.set_acc(this->mepsrel.imag(), this->mepsabs.imag());
            try {
                mquad.integrate_sing(im_func, a, b, im, err_im, limit);
            } catch(std::runtime_error &e) {
                res = re + cx_double(0.,1.) * im;
                abserr = err_re + cx_double(0.,1.) * err_im;
                throw;
            }
            throw;
        }
        mquad.set_acc(this->mepsrel.imag(), this->mepsabs.imag());
        try {
            mquad.integrate_sing(im_func, a, b, im, err_im, limit);
        } catch(std::runtime_error &e) {
            res = re + cx_double(0.,1.) * im;
            abserr = err_re + cx_double(0.,1.) * err_im;
            throw;
        }
        res = re + cx_double(0.,1.) * im;
        abserr = err_re + cx_double(0.,1.) * err_im;
    }

    cx_double
    integrate_sing(const Func& func, double a, double b, cx_double& abserr, std::size_t limit = 2000) {
        cx_double res;
        integrate_sing(func, a, b, res, abserr, limit);
        return res;
    }


    using base::set_acc;
    void set_acc(double epsrel, double epsabs) { base::set_acc(cx_double(epsrel, epsrel), cx_double(epsabs, epsabs)); }

protected:
    template<typename Tf>
    struct Wrapper_t {
        explicit Wrapper_t(const Tf& f): mpf(&f) { }
        virtual ~Wrapper_t() = default;
        virtual double operator()(double x) const = 0;
        const Tf *mpf;
    };
    template<typename Tf>
    struct Re_t: Wrapper_t<Tf> {
        explicit Re_t(const Tf& f): Wrapper_t<Tf>(f) { }
        virtual double operator()(double x) const override {
            return std::real( (*this->mpf)(x) );
        }
    };
    template<typename Tf>
    struct Im_t: Wrapper_t<Tf> {
        explicit Im_t(const Tf& f): Wrapper_t<Tf>(f) { }
        virtual double operator()(double x) const override {
            return std::imag( (*this->mpf)(x) );
        }
    };

    Rag_gsl<Wrapper_t<Func>> mquad;

    template<typename Ta, typename Tb>
    int integrate_impl(Ta a, Tb b, cx_double& res, cx_double& abserr, std::size_t limit) { 
        Re_t<Func> re_func(*this->mpfunc);
        Im_t<Func> im_func(*this->mpfunc);
        double re, im, err_re, err_im;

        mquad.set_acc(this->mepsrel.real(), this->mepsabs.real());
        try {
            mquad.integrate(re_func, a, b, re, err_re, limit);
        } catch (std::runtime_error& e) {
            mquad.set_acc(this->mepsrel.imag(), this->mepsabs.imag());
            try {
                mquad.integrate(im_func, a, b, im, err_im, limit);
            } catch(std::runtime_error &e) {
                res = re + cx_double(0.,1.) * im;
                abserr = err_re + cx_double(0.,1.) * err_im;
                throw;
            }
            throw;
        }
        mquad.set_acc(this->mepsrel.imag(), this->mepsabs.imag());
        try {
            mquad.integrate(im_func, a, b, im, err_im, limit);
        } catch(std::runtime_error &e) {
            res = re + cx_double(0.,1.) * im;
            abserr = err_re + cx_double(0.,1.) * err_im;
            throw;
        }
        res = re + cx_double(0.,1.) * im;
        abserr = err_re + cx_double(0.,1.) * err_im;
        return 0;
    }

    virtual int
    integrate_impl(double a, std::true_type, double b, std::true_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return integrate_impl(a, b, res, abserr, limit);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::false_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return integrate_impl(mquad.neg_inf(), mquad.pos_inf(), res, abserr, limit);
    }
    virtual int
    integrate_impl(double a, std::true_type, double b, std::false_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return integrate_impl(a, mquad.pos_inf(), res, abserr, limit);
    }
    virtual int
    integrate_impl(double a, std::false_type, double b, std::true_type,
                   cx_double& res, cx_double& abserr, std::size_t limit) override {
        return integrate_impl(mquad.neg_inf(), b, res, abserr, limit);
    }
};

//**************************************************************************************

template<typename Res, typename Func>
class Bone_glfixed: public Bone<Res, Func> {
    using base = Bone<Res, Func>;
public:
    Bone_glfixed() = delete;
    explicit Bone_glfixed(std::size_t nPoints):
        base(), mnPoints(nPoints) { }

    Res integrate(const Func& func, double a, double b) {
        this->mpfunc = &func;
        return integrate_impl(a, b);
    }

protected:
    std::size_t mnPoints;

    static Res
    func_wrapper(double x, void *pq) {
        auto pquad = static_cast<Bone_glfixed *>(pq);
        return (*pquad->mpfunc)(x);
    }

    virtual Res
    integrate_impl(double a, double b) { return 0.; };
};

template<typename Func = std::complex<double>(double)>
class Cglfixed: public Bone_glfixed<std::complex<double>, Func> {
protected:
    using cx_double = std::complex<double>;
    using base = Bone_glfixed<cx_double, Func>;
public:
    explicit Cglfixed(std::size_t nPoints = 128):
        base(nPoints), mtable(glfixed_table_alloc(nPoints), glfixed_table_free)
    {
        if (mtable.get() == nullptr) {
            gsl_set_error_handler(this->merr_h);
            throw std::bad_alloc{};
        }
        mf.function = &this->func_wrapper;
        mf.params = this;
        gsl_set_error_handler(this->merr_h);
    }

    Cglfixed(const Cglfixed&) = delete;
    Cglfixed(Cglfixed&&) = default;
    Cglfixed& operator=(const Cglfixed&) = delete;
    Cglfixed& operator=(Cglfixed&&) = default;
    
    using base::integrate;

    virtual void
    integrate(const Func& func, double a, double b, cx_double& res, cx_double& untouched, std::size_t nPoints = 128)
    override {
        if (nPoints != this->mnPoints) mtable.reset(glfixed_table_alloc(nPoints));
        res = base::integrate(func, a, b);
    }

    virtual cx_double
    integrate(const Func& func, double a, double b, cx_double& untouched, std::size_t nPoints = 128) override {
        cx_double res;
        integrate(func, a, b, res, untouched, nPoints);
        return res;
    }

protected:
    using table_t = std::unique_ptr<glfixed_table, decltype(&glfixed_table_free)>;
    table_t mtable;
    cr_function mf;

    virtual cx_double
    integrate_impl(double a, double b) override {
        return cr_gsl_integration_glfixed(&this->mf, a, b, mtable.get());
    }
};

template<typename Func = double(double)>
class Rglfixed: public Bone_glfixed<double, Func> {
protected:
    using base = Bone_glfixed<double, Func>;
public:
    explicit Rglfixed(std::size_t nPoints = 128):
        base(nPoints), mtable(gsl_integration_glfixed_table_alloc(nPoints), gsl_integration_glfixed_table_free)
    {
        if (mtable.get() == nullptr) {
            gsl_set_error_handler(this->merr_h);
            throw std::bad_alloc{};
        }
        mf.function = &this->func_wrapper;
        mf.params = this;
        gsl_set_error_handler(this->merr_h);
    }

    Rglfixed(const Rglfixed&) = delete;
    Rglfixed(Rglfixed&&) = default;
    Rglfixed& operator=(const Rglfixed&) = delete;
    Rglfixed& operator=(Rglfixed&&) = default;
    
    using base::integrate;

    virtual void
    integrate(const Func& func, double a, double b, double& res, double& untouched, std::size_t nPoints = 128)
    override {
        if (nPoints != this->mnPoints) mtable.reset(gsl_integration_glfixed_table_alloc(nPoints));
        return base::integrate(func, a, b);
    }

    virtual double
    integrate(const Func& func, double a, double b, double& untouched, std::size_t nPoints = 128) override {
        double res;
        integrate(func, a, b, res, untouched, nPoints);
        return res;
    }

protected:
    using table_t = std::unique_ptr<gsl_integration_glfixed_table, decltype(&gsl_integration_glfixed_table_free)>;
    table_t mtable;
    gsl_function mf;

    virtual double
    integrate_impl(double a, double b) override {
        return gsl_integration_glfixed(&this->mf, a, b, mtable.get());
    }
};

//**************************************************************************************

} //namespace quad_1d

#endif
