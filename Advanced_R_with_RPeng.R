library(swirl)
library(purrr)                 

# MAP() map_lgl(), map_chr(), and map_dbl() functions

map_chr(c(5, 4, 3, 2, 1), function(x){
  c("one", "two", "three", "four", "five")[x]
})
 

map_lgl(c(1, 2, 3, 4, 5), function(x){
  x > 3
})


map_if(1:5, function(x){
    x %% 2 == 0
  },
  function(y){
    y^2
  }) %>% unlist()


#The map_at() function only applies the provided function to elements of a vector specified by their indexes. map_at() always returns a list 
map_at(seq(100, 500, 100), c(1, 3, 5), function(x){
  x - 10
}) %>% unlist()

#The first two arguments should be two vectors of the same length, followed by a function which will be evaluated with an element of the first vector as the first argument and an element of the second vector as the second argument.
map2_chr(letters, 1:26, paste)

#The pmap() family of functions is similar to map2(), however instead of mapping across two vectors or lists, you can map across any number of lists. 
pmap_chr(list(
  list(1, 2, 3),
  list("one", "two", "three"),
  list("uno", "dos", "tres")
), paste)



#REDUCE()

#List or vector reduction iteratively combines the first element of a vector with the second element of a vector, then that combined result is combined with the third element of the vector, and so on until the end of the vector is reached. The function to be applied should take at least two arguments.



reduce(c(1, 3, 5, 7), function(x, y){
  message("x is ", x)
  message("y is ", y)
  message("")
  x + y
})


reduce(letters[1:4], function(x, y){
  message("x is ", x)
  message("y is ", y)
  message("")
  paste0(x, y)
}, .dir = "backward")


#SEARCH
library(tidyselect)

"a" %in% letters


#The detect() function takes a vector and a predicate function as arguments and it returns the first element of the vector for which the predicate function returns TRUE:
detect(20:40, function(x){
  x > 22 && x %% 2 == 0
})

detect_index(20:40, function(x){
  x > 22 && x %% 2 == 0
})


#FILTER

keep(1:20, function(x){
  x %% 2 == 0
})

#opposite, it leaves only those one. who doesn't satisfy
discard(1:20, function(x){
  x %% 2 == 0
})

#The every() function returns TRUE only if every element in the vector satisfies the predicate function, 
every(1:20, function(x){
  x %% 2 == 0
})

#while the some() function returns TRUE if at least one element in the vector satisfies the predicate function:
some(1:20, function(x){
  x %% 2 == 0
})


#COMPOSE

#the compose() function combines any number of functions into one function:
n_unique <- compose(length, unique)
# The composition above is the same as:
# n_unique <- function(x){
#   length(unique(x))
# }

rep(1:5, 1:5)

n_unique(rep(1:5, 1:5))


#PARTIAL APPLICATION

library(purrr)

mult_three_n <- function(x, y, z){
  x * y * z
}

#Partial application of functions can allow functions to behave a little like data structures. Using the partial() function from the purrrpackage you can specify some of the arguments of a function, and then partial() will return a function that only takes the unspecified arguments. 

mult_by_15 <- partial(mult_three_n, x = 3, y = 5)

mult_by_15(z = 4)

#By using partial application you can bind some data to the arguments of a function before using that function elsewhere.


#SIDE EFFECTS AND WALK

# If you want to evaluate a function across a data structure you should use the walk() function from purrr.

walk(c("Friends, Romans, countrymen,","lend me your ears;","I come to bury Caesar,", 
       "not to praise him."), message)

#RECURSION

vector_sum_rec <- function(v){
  if(length(v) == 1){
    v
  } else {
    v[1] + vector_sum_rec(v[-1])
  }
}

#Let’s write a function to compute the nth digit of the Fibonacci sequence such that the first number in the sequence is 0, the second number is 1, and then all proceeding numbers are the sum of the n - 1 and the n - 2 Fibonacci number. 

fib <- function(n){
  stopifnot(n > 0)
  if(n == 1){
    0
  } else if(n == 2){
    1
  } else {
    fib(n - 1) + fib(n - 2)
  }
}

map_dbl(1:12, fib)



fib_tbl <- c(0, 1, rep(NA, 23))

fib_mem <- function(n){
  stopifnot(n > 0)
  
  if(!is.na(fib_tbl[n])){
    fib_tbl[n]
  } else {
    fib_tbl[n - 1] <<- fib_mem(n - 1)
    fib_tbl[n - 2] <<- fib_mem(n - 2)
    fib_tbl[n - 1] + fib_tbl[n - 2]
  }
}

map_dbl(1:12, fib_mem)






library(purrr)
library(microbenchmark)
library(tidyr)
library(magrittr)
library(dplyr)

fib_data <- map(1:10, function(x){microbenchmark(fib(x), times = 100)$time})
names(fib_data) <- paste0(letters[1:10], 1:10)
fib_data <- as.data.frame(fib_data)

fib_data %<>%
  gather(num, time) %>%
  group_by(num) %>%
  summarise(med_time = median(time))

memo_data <- map(1:10, function(x){microbenchmark(fib_mem(x))$time})
names(memo_data) <- paste0(letters[1:10], 1:10)
memo_data <- as.data.frame(memo_data)

memo_data %<>%
  gather(num, time) %>%
  group_by(num) %>%
  summarise(med_time = median(time))

plot(1:10, fib_data$med_time, xlab = "Fibonacci Number", ylab = "Median Time (Nanoseconds)",
     pch = 18, bty = "n", xaxt = "n", yaxt = "n")
axis(1, at = 1:10)
axis(2, at = seq(0, 350000, by = 50000))
points(1:10 + .1, memo_data$med_time, col = "blue", pch = 18)
legend(1, 300000, c("Not Memorized", "Memoized"), pch = 18, 
       col = c("black", "blue"), bty = "n", cex = 1, y.intersp = 1.5)


# EXPRESSIONS

two_plus_two <- quote(2+2)
eval(two_plus_two)

string_example <- "2 + 2"
string_to_expression <- parse( text = string_example)
eval(string_to_expression)

# opposite of quote is deparse()

deparse(two_plus_two)

# How to modify the functions to be executed

sum_expression <- quote(sum(1,5))
eval(sum_expression)
sum_expression[[1]]

sum_expression[[1]] <- quote(paste0)

eval(sum_expression)


# CALL() You can compose expressions using the call() function. T
# he first argument is a string containing the name of a function, followed by the arguments 
# that will be provided to that function.

sum_expression_2 <- call("sum", 40,50)
eval(sum_expression_2)


return_expression <- function(...){
  match.call()
}

return_expression(2, col = "blue", FALSE)

first_arg <- function(...){
  expr <- match.call()
  first_arg_expr <- expr[[2]]
  first_arg <- eval(first_arg_expr)
  if(is.numeric(first_arg)){
    paste("The first argument is", first_arg)
  } else {
    "The first argument is not numeric."
  }
}

first_arg(2, 4, "seven", FALSE)
first_arg("two", 4, "seven", FALSE)


# ENVIRONMENTS

my_new_env <- new.env()
my_new_env$x <- 4
my_new_env$x

assign("y", 9, envir = my_new_env)
get("y", envir = my_new_env)

my_new_env$y

ls(my_new_env)

rm(y, envir = my_new_env)
exists("y", envir = my_new_env)
exists("x", envir = my_new_env)

my_new_env$x
my_new_env$y

search()

# EXECUTION

x <- 10

my_func <- function(){
  x <- 5
  return(x)
}

my_func()


another_func <- function(){
  return(x)
}

another_func()

# complex assignment operator <<- , 
# also use <<- to assign names to values that have not been yet been defined in the global environment from inside a function:

x <- 10
x

assign1 <- function(){
  x <<- "Wow!"
}

assign1()
x

a_variable_name
exists("a_variable_name")

assign2 <- function(){
  a_variable_name <<- "Magic!"
}

assign2()
exists("a_variable_name")
a_variable_name


# Generating Errors

stop("Something erroneous has occured!")
name_of_function <- function(){
  stop("Something bad happened.")
}

name_of_function()

# stopifnot()

error_if_n_is_greater_than_zero <- function(n){
  stopifnot(n <= 0)
  n
}

error_if_n_is_greater_than_zero(5)

# warning()
warning("This is a warning message")
make_NA <- function(x){
  warning("Generating an NA.")
  NA
}

make_NA()

#HANDLING ERRORS

beera <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
           },
           warning = function(w){
             message("A warning occured:\n", w)
           },
           finally = {
             message("Finally done!")
           })
}


beera(2+2)


### FUNCTIONAL PROGRAMMING TASK FOR 3RD WEEK
int_to_string(4)
gt(a, b)
is_even(x)
square(x)
add_talk(x, y)
paste_talk(x, y)


### DEBUGGING


# traceback()
check_n_value <- function(n) {
  if(n > 0) {
    stop("n should be <= 0")
  }
}
error_if_n_is_greater_than_zero <- function(n){
  check_n_value(n)
  n
}
error_if_n_is_greater_than_zero(5)
traceback()


# browser()
check_n_value <- function(n) {
  if(n > 0) {
    browser()  ## Error occurs around here
    stop("n should be <= 0")
  }
}
error_if_n_is_greater_than_zero(5)

# trace() function to make temporary code modifications.
trace("check_n_value")

error_if_n_is_greater_than_zero(5)

as.list(body(check_n_value))

as.list(body(check_n_value)[[2]])

as.list(body(check_n_value)[[2]])[3]

# Now we can see the call to stop() is the third sub-expression within the second 
# expression of the overall function. We can specify this to trace() by passing an integer 
# vector wrapped in a list to the at argument.
trace("check_n_value", browser, at = list(c(2, 3)))

# The trace() function has a side effect of modifying the function and converting into 
# a new object of class “functionWithTrace”.

check_n_value

body(check_n_value)

trace("check_n_value", quote({
  if(n == 5) {
    message("invoking the browser")
    browser()
  }
}), at = 2)
body(check_n_value)


# debugging functions within a package is another key use case for trace(). 

# if we wanted to insert tracing code into the glm() function within the stats package, 
# the only addition to the trace() call we would need is to provide the namespace 
# information via the where argument.

trace("glm", browser, at = 4, where = asNamespace("stats"))
body(stats::glm)[1:5]



# Using debug() and debugonce()























