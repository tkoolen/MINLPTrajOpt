module BARONTools

export
    parse_sum

using DataFrames

@enum SubProblemType begin
    RELAXATION
    PROBING
end

const SUM_FILE_DATA_HEADER = "  Iteration    Open nodes         Time (s)    Lower bound      Upper bound\n"

macro dataframe(args...)
    argexprs = map(arg -> Expr(:kw, arg, :($(esc(arg)))), args)
    :(DataFrame($(argexprs...)))
end

function sum_iteration_data_lines(io::IO)
    mark(io)
    N = 0
    lastchar = ' '
    while !eof(io)
        c = read(io, Char)
        if c == '\n'
            if lastchar == '\n'
                break
            else
                N += 1
            end
        end
        lastchar = c
    end
    reset(io)
    N
end

function readuntil_pred(pred, io::IO)
    buf = IOBuffer()
    while !eof(io)
        c = read(io, Char)
        if pred(c)
            nbytes = write(devnull, c)
            skip(io, -nbytes)
            break
        end
        write(buf, c)
    end
    return String(take!(buf))
end

function parse_sum(io::IO)
    header = readuntil(io, SUM_FILE_DATA_HEADER)
    N = sum_iteration_data_lines(io)
    iteration = Vector{Int}(undef, N)
    subproblem_type = Vector{SubProblemType}(undef, N)
    open_nodes = Vector{Int}(undef, N)
    time = Vector{Float64}(undef, N)
    lower_bound = Vector{Float64}(undef, N)
    upper_bound = Vector{Float64}(undef, N)

    for i in 1 : N
        # iteration
        skipchars(x -> x === '*', io)
        skipchars(isspace, io)
        iteration[i] = parse(Int, readuntil_pred(!isnumeric, io))

        # subproblem_type
        subproblem_type[i] = read(io, Char) == '*' ? PROBING : RELAXATION

        # open_nodes
        skipchars(isspace, io)
        open_nodes[i] = parse(Int, readuntil_pred(!isnumeric, io))

        # time
        skipchars(isspace, io)
        time[i] = parse(Float64, readuntil_pred(isspace, io))

        # lower_bound
        skipchars(isspace, io)
        lower_bound[i] = parse(Float64, readuntil_pred(isspace, io))

        # upper_bound
        skipchars(isspace, io)
        upper_bound[i] = parse(Float64, readuntil_pred(isspace, io))
    end
    @dataframe iteration subproblem_type open_nodes time lower_bound upper_bound
end

parse_sum(sumfile::AbstractString) = open(parse_sum, sumfile)

end # module
