#include "errorvisitor.hpp"

void ErrorVisitor::push_error(Expression *e) {
    if (e->has_error()) {
        has_error_ = true;

        if (record_) {
            record_->errors().push_back({e->error_message(), e->location()});
        }

        if (!quiet_) {
            std::cout
                << red("error: ") << white(pprintf("%", e->location()))
                << "\n    " << e->error_message() << "\n";

        }
    }

    if (e->has_warning()) {
        if (record_) {
            record_->warnings().push_back({e->warning_message(), e->location()});
        }

        if (!quiet_) {
            std::cout
                << purple("warning: ") << white(pprintf("%", e->location()))
                << "\n    " << e->warning_message() << "\n";
        }
    }
}

/*
 * we use a post order walk to print the erros in an expression after those
 * in all of its children
 */

void ErrorVisitor::visit(Expression *e) {
    push_error(e);
}

// traverse the statements in a procedure
void ErrorVisitor::visit(ProcedureExpression *e) {
    for(auto& expression : e->args()) {
        expression->accept(this);
    }

    e->body()->accept(this);
    push_error(e);
}

// traverse the statements in a function
void ErrorVisitor::visit(FunctionExpression *e) {
    for(auto& expression : e->args()) {
        expression->accept(this);
    }

    e->body()->accept(this);
    push_error(e);
}

// an if statement
void ErrorVisitor::visit(IfExpression *e) {
    e->true_branch()->accept(this);
    if(e->false_branch()) {
        e->false_branch()->accept(this);
    }

    push_error(e);
}

void ErrorVisitor::visit(BlockExpression* e) {
    for(auto& expression : e->statements()) {
        expression->accept(this);
    }

    push_error(e);
}

void ErrorVisitor::visit(InitialBlock* e) {
    for(auto& expression : e->statements()) {
        expression->accept(this);
    }

    push_error(e);
}

// unary expresssion
void ErrorVisitor::visit(UnaryExpression *e) {
    e->expression()->accept(this);
    push_error(e);
}

// binary expresssion
void ErrorVisitor::visit(BinaryExpression *e) {
    e->lhs()->accept(this);
    e->rhs()->accept(this);
    push_error(e);
}

// binary expresssion
void ErrorVisitor::visit(CallExpression *e) {
    for(auto& expression: e->args()) {
        expression->accept(this);
    }
    push_error(e);
}

