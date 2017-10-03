#pragma once

namespace arbor
{

class XArbor : public std::exception
    {
    public:
                            XArbor() throw() {}
                            XArbor(const std::string s) throw() : _msg() {_msg = s;}
        virtual             ~XArbor() throw() {}
        const char *        what() const throw() {return _msg.c_str();}
        
    private:
    
        std::string         _msg;
    };
    
}