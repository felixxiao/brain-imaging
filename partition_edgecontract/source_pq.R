MinPriorityQueue = R6Class('MinPriorityQueue',
  public = list(
    priority = numeric(0),
    item = NA,
    N = 0,

    is_empty = function() { self$N == 0 },

    exchange = function(i, j)
    {
      tmp = self$priority[i]
      self$priority[i] = self$priority[j]
      self$priority[j] = tmp

      tmp = self$item[i]
      self$item[i] = self$item[j]
      self$item[j] = tmp
    },

    swim = function(k)
    {
      while (k > 1)
      {
        if (self$priority[k%/%2] < self$priority[k]) break
        self$exchange(k, k %/% 2)
        k = k%/%2
      }
    },

    sink = function(k)
    {
      while (2*k <= self$N)
      {
        j = 2 * k +
            (2*k < self$N & self$priority[2*k + 1] < self$priority[2*k])
        if (self$priority[k] < self$priority[j]) break

        self$exchange(k, j)
        k = j
      }
    },

    insert = function(item, priority)
    {
      self$N = self$N + 1
      self$priority[self$N] = priority
      self$item[self$N] = item
      self$swim(self$N)
    },

    remove_min = function()
    {
      if (self$N == 0)
        stop('Empty queue')
      item = self$item[1]
      priority = self$priority[1]
      self$exchange(1, self$N)
      self$priority[self$N] = NA
      self$item[self$N] = NA
      self$N = self$N - 1
      self$sink(1)
      return(list(item = item, priority = priority))
    }
  )
)

"
pq = MinPriorityQueue$new()
pq$insert('d', 4)
pq$insert('b', 2)
pq$insert('e', 5)
pq$insert('c', 3)
pq$insert('a', 1)
pq$insert('f', 6)

while (! pq$is_empty())
  cat(pq$remove_min()$item, ' ')
"