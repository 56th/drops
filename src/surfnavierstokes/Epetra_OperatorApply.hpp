//
// Created by Alexander Zhiliakov alex@math.uh.edu on 10/18/19.
//

#ifndef EPETRA_OPERATOR_APPLY_HPP
#define EPETRA_OPERATOR_APPLY_HPP

#include <functional>
#include <string>
#include "Epetra_Operator.h"

namespace DROPS {

    class Epetra_OperatorApply : public Epetra_Operator {
    public:
        using ApplyType = std::function<void(Epetra_MultiVector const &, Epetra_MultiVector&)>;
    private:
        ApplyType apply;
    public:
        Epetra_OperatorApply(ApplyType const & aplly) : apply(apply) {}
        virtual ~Epetra_OperatorApply() {}
        int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const final {
            apply(X, Y);
            return 0;    
        }
        int SetUseTranspose(bool UseTranspose) final {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        double NormInf() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        const char *Label() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        bool UseTranspose() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        bool HasNormInf() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        const Epetra_Comm &Comm() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        const Epetra_Map &OperatorDomainMap() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
        const Epetra_Map &OperatorRangeMap() const {
            std::string funcName = __func__;
            throw std::logic_error(funcName + ": not implemented");
        }
    };
}

#endif //DROPS_BELOSDROPSADAPTER_HPP
