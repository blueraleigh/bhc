#library(bhc)
#data(snakediet, package='macroevolution')

#x = xtabs(snakediet$count ~ snakediet$species + snakediet$food)


#res = bhc.multinomial(x)


#snakediet = read.table("/Users/mgrundler/Dropbox/manuscripts/snakediet_ms/data/snakediet.txt",
#    header=TRUE)


#x = xtabs(snakediet$count ~ snakediet$predator + snakediet$prey)


#system.time(res <- bhc.multinomial(x))

#g = cutree(res, 53)

#state = 46
#barplot(colSums(x[names(which(g==state)), ]) / sum(x[names(which(g==state)), ]),
#    cex.names=0.35)



#foo = function(node, res, alpha=0.5) {
#    tips = function(node) {
#        ans = c()
#        foo = function(node) {
#            if (node > 0) {
#                lf = res$merge[node, 1]
#                rt = res$merge[node, 2]
#                foo(lf)
#                foo(rt)
#            } else {
#                ans <<- c(ans, node)
#            }
#        }
#        foo(node)
#        -ans
#    }

#    bar = function(node) {
#        if (node > 0) {
#            if (exp(res$height[node]) < alpha) {
#                lf = res$merge[node, 1]
#                rt = res$merge[node, 2]
#                bar(lf)
#                bar(rt)
#            } else {
#                g <<- g + 1
#                #ans[[g]] <<- tips(node)
#                ans[tips(node)] <<- g
#            }
#        } else {
#            g <<- g + 1
#            #ans[[g]] <<- -node
#            ans[tips(node)] <<- g
#        }
#    }
#    #ans = list()
#    ans = integer(node+1)
#    g = 0
#    bar(node)
#    ans
#}


#g = foo(881, res, .5)


#rownames(x)[tips(804, res)]
