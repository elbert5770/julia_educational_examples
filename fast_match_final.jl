using BenchmarkTools

function sort_with_tolerance(s1,rows,eps)
    for i in range(1,rows-1)
        if abs(s1[i,1] - s1[i+1,1]) < eps
            s1[i,3] = s1[i+1,1]
        else
            if abs(s1[i,2] - s1[i+1,2]) < eps
                s1[i,2] = s1[i+1,2]
            end
        end
    end
end

function match(cell_center1,cell_center2,rows1,rows2,cell_match1,cell_match2,cell1,cell2,eps)
    # Sort ararys by values in col3, col2, then col1
    rcell_center1 = reverse(cell_center1,dims=2)
    sort_with_tolerance(rcell_center1,rows1,eps)
    p1 = sortperm(eachslice(rcell_center1; dims=1))
    rcell_center2 = reverse(cell_center2,dims=2)  
    sort_with_tolerance(rcell_center2,rows2,eps)
    p2 = sortperm(eachslice(rcell_center2; dims=1))

    # Create sorted arrays s1 & s2
    s1 = cell_center1[p1,:]
    s2 = cell_center2[p2,:]
    # @show s1
    # @show s2
    m = rows1
    n = rows2
    i = 1
    j = 1
    
    count = 0
    counter = 0
    stemp = Vector{Float64}(undef,3)
    stemp_abs = Vector{Int64}(undef,3)
    stemp_sign = 0
    while i <= m && j <= n
        # Do all three coordinates match? Record in cell_match else advance array with coordinate closest to origin
        if abs(s1[i,3] - s2[j,3])  < eps && abs(s1[i,2] - s2[j,2])  < eps && abs(s1[i,1] - s2[j,1])  < eps
            # One of the two matches will be the owner for distributed processing later
            owner = rand((0, 1))

            ## 656 ms, 16 allocations, 91.55 MiB
            cell_match1[p1[i],1:4] .= owner,cell2,p2[j],p1[i]
            cell_match2[p2[j],1:4] .= 1-owner,cell1,p1[i],p2[j]
            # @show i,j,s1[i,:],s2[j,:]
            ## 642 ms, 7690 allocations, 92.26 MiB
            # cell_match1[p1[i],1:4] = [owner,cell2,p2[j],p1[i]]
            # cell_match2[p2[j],1:4] = [1-owner,cell1,p1[i],p2[j]]

            ## 644 ms, 7666 allocations, 92 MiB
            # cell_match1[p1[i],1:4] .= [owner,cell2,p2[j],p1[i]]
            # cell_match2[p2[j],1:4] .= [1-owner,cell1,p1[i],p2[j]]

            ## 629 ms, 7884 allocations, 92.27 MiB
            # @views cell_match1[p1[i],1:4] = [owner,cell2,p2[j],p1[i]]
            # @views cell_match2[p2[j],1:4] = [1-owner,cell1,p1[i],p2[j]]

            ## 640 ms, 16 allocations, 91.55 MiB
            # cell_match1[p1[i],1] = owner
            # cell_match1[p1[i],2] = cell2
            # cell_match1[p1[i],3] = p2[j]
            # cell_match1[p1[i],4] = p1[i]
            # cell_match2[p2[j],1] = 1-owner
            # cell_match2[p2[j],2] = cell1
            # cell_match2[p2[j],3] = p1[i]
            # cell_match2[p2[j],4] = p2[j]
            i += 1
            j += 1
        else
            # Need to identify which is the first difference between s1 and s2 in the order
            # of z, y, x.  If the value of the element is greater in s1, then increment j
            # else increment i

            ## 664.751 ms (18 allocations: 91.55 MiB)
            @views stemp .= (s1[i,1:3] .- s2[j,1:3])
            stemp_abs .= abs.(stemp) .> eps # 0 if false, so it is within eps of zero
            count = stemp_abs[3] > 0 ? 3 : (stemp_abs[2] > 0 ? 2 : 1)            
            stemp_sign = sign(stemp[count])/2
            # @show i,j,s1[i,:],s2[j,:],stemp,stemp_abs,stemp_sign,count,stemp_sign,sign(stemp[count])
            i += round(Int64,(1-stemp_sign)/2)
            j += round(Int64,(stemp_sign+1)/2)
          
                
            # 666.481 ms (16 allocations: 91.55 MiB)
            # if s1[i,3] - s2[j,3] < -eps
            #     i += 1
            # else
            #     if s1[i,3] - s2[j,3] > eps
            #         j += 1
            #     else
            #         if s1[i,2] - s2[j,2] < -eps
            #             i += 1
            #         else
            #             if s1[i,2] - s2[j,2] > eps
            #                 j += 1
            #             else
            #                 if s1[i,1] - s2[j,1] < -eps
            #                     i += 1
            #                 else
            #                     if s1[i,1] - s2[j,1] > eps
            #                         j += 1
            #                     end
            #                 end
            #             end
            #         end

            #     end
            # end
        end
        
    
    end
    return nothing
end

function main()
    eps = 1e-6
    cell1 = 123
    cell2 = 254
    rows1 = 1000000
    rows2 = 500000
    cols = 3
    cell_center1 = rand(0.0:500.0,(rows1,cols)) .+ rand(-1e-8:1e-8,(rows1,cols))
    cell_center2 = rand(0.0:500.0,(rows2,cols)) .+ rand(-1e-8:1e-8,(rows2,cols))
    row1 = 1
    row2 = 1
    num_forced_matches = 0
    forced_matches = Array{Int64,2}(undef,num_forced_matches,2)

    # If arrays are small, no matches will be found.  Ensure that at least 15 matches are found
    for i in range(1,num_forced_matches)
        row1 = rand(1:rows1)
        row2 = rand(1:rows2)
        forced_matches[i,1:2] .= row1,row2
        cell_center2[row2,:] .= cell_center1[row1,:]
        # @show row1,row2,cell_center1[row1,:],cell_center2[row2,:]
    end

    # Allocate arrays for matches
    cell_match1 = Array{Int64,2}(undef,rows1,4)*0
    cell_match2 = Array{Int64,2}(undef,rows2,4)*0

    match(cell_center1,cell_center2,rows1,rows2,cell_match1,cell_match2,cell1,cell2,eps)
    @btime match($cell_center1,$cell_center2,$rows1,$rows2,$cell_match1,$cell_match2,$cell1,$cell2,$eps)

    # @show sortslices(cell_match2,dims=1)

end

main()