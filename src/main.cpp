#include "trigonometry.hpp"
#include <iostream>

int main()
{
    std::string answer;
    while (answer != std::to_string(3))
    {
        std::cout << "Enter \n1 - to calculate sin \n2 - to calculate cos \n3 - to exit" << std::endl;
        std::cin >> answer;
        if (answer == std::to_string(1))
        {
            std::string arg;
            std::cout << "Enter argument in degrees:";
            std::cin >> arg;
            std::cout << "Answer: " << Trigonometrix::sinDeg<double,true>(atof(arg.c_str())) << std::endl;
        }
        else if (answer == std::to_string(2))
        {
            std::string arg;
            std::cout << "Enter argument in degrees:";
            std::cin >> arg;
            std::cout << "Answer: " << Trigonometrix::cosDeg<double,true>(atof(arg.c_str())) << std::endl;
        }
    }
    return 0;
}
