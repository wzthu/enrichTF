#' @name EnrichStep-class
#' @rdname  EnrichStep-class
#' @title Base class of this package
#' @description This class is inherit from \code{Step} in pipeFrame package,
#' no more method is extended or override. Please see \code{Step} class for detail.
#' @export
setClass(Class = "EnrichStep",
         contains = "Step"
)


#' @importFrom pipeFrame input "input<-" output "output<-" "param" "param<-",
