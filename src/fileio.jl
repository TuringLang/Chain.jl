#################### File I/O ####################

function readcoda(output::AbstractString, index::AbstractString)
  out = readdlm(output, Any)
  ind = readdlm(index, Any)

  firstind = ind[:, 2]
  firstiter = out[firstind, 1]
  lastind = ind[:, 3]
  lastiter = out[lastind, 1]

  thin = Int((lastiter[1] - firstiter[1]) / (lastind[1] - firstind[1]))
  window = maximum(firstiter):thin:minimum(lastiter)
  startind = firstind + (first(window) - firstiter) / step(window)
  stopind = lastind - (lastiter - last(window)) / step(window)

  names = AbstractString[ind[:, 1]...]

  value = Array{Float64}(length(window), length(names))
  for i in 1:size(value, 2)
    inds = Int(startind[i]):Int(stopind[i])
    value[:, i] = out[inds, 2]
  end

  Chains(value; iterations=window, names=names)
end
