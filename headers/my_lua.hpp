#ifndef MY_LUA_HPP
#define MY_LUA_HPP

#include <iostream>

namespace lua
{
extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

struct lua_value_wrapper
{
   lua_State* L ;

   lua_value_wrapper(std::string const& script_file) :
      L(lua_open()) 
   {
      // load the libs
      luaL_openlibs(L)  ;
 
      //run a Lua scrip here
      luaL_dofile(L, script_file.c_str()) ;
   }
   
   ~lua_value_wrapper()
   {
      lua_close(L);
   }
};

template<typename TargetType>
inline TargetType extract(lua_value_wrapper& lua
                          , std::string const& variable_name)
{}

template<>
inline double extract<double>(lua_value_wrapper& lua
                              , std::string const& variable_name)
{
   lua_getglobal(lua.L, variable_name.c_str());
   if (!lua_isnumber(lua.L, -1)) {
      std::cerr << "erroneous conversion to double for: "
                << variable_name 
                << std::endl ;
      exit (-1) ;
      return 0.;
   }
   else
   {
      return lua_tonumber(lua.L, -1) ;
   }
}

template<>
inline long extract<long>(lua_value_wrapper& lua
                          , std::string const& variable_name)
{
   lua_getglobal(lua.L, variable_name.c_str());
   if (!lua_isnumber(lua.L, -1)) {
      std::cerr << "erroneous conversion to long for: "
                << variable_name
                << std::endl ;
      exit (-1) ;
      return 0;
   }
   else
   {
      return lua_tointeger(lua.L, -1) ;
   }
}

template<>
inline std::string extract<std::string>(lua_value_wrapper& lua
                                        , std::string const& variable_name)
{
   lua_getglobal(lua.L, variable_name.c_str());
   if (!lua_isstring(lua.L, -1)) {
      std::cerr << "erroneous conversion to string for: "
                << variable_name
                << std::endl ;
      exit (-1) ;
      return std::string();
   }
   else
   {
      return std::string(lua_tostring(lua.L, -1)) ;
   }


   return 0;
}



}


#endif // MY_LUA_HPP
