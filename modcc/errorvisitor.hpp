#pragma once

#include <iostream>
#include "visitor.hpp"
#include "expression.hpp"

class ErrorVisitor : public Visitor {
public:
    ErrorVisitor(error_stack* record, bool quiet = true)
        : record_(record), quiet_(quiet)
    {}

    void visit(Expression *e)           override;
    void visit(ProcedureExpression *e)  override;
    void visit(FunctionExpression *e)   override;
    void visit(UnaryExpression *e)      override;
    void visit(BinaryExpression *e)     override;
    void visit(CallExpression *e)       override;

    void visit(BlockExpression *e)      override;
    void visit(InitialBlock *e)         override;
    void visit(IfExpression *e)         override;

    operator bool() const { return has_error_; }

private:
    void push_error(Expression*);

    error_stack* record_ = nullptr;
    bool quiet_ = true;
    bool has_error_ = false;
};

inline bool collect_errors(Expression* e, error_stack* record = nullptr, bool quiet = true) {
    ErrorVisitor v(record, quiet);
    v.visit(e);
    return v;
}

inline bool collect_errors(expression_ptr& e, error_stack* record = nullptr, bool quiet = true) {
    return collect_errors(e.get(), record, quiet);
}
