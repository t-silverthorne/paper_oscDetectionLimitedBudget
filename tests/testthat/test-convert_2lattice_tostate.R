test_that("compare to ref", {
lat1=c(1:4)/4-1/4
lat2=c(1:4)/4-1/4

lat_ref = c(1:8)/8-1/8
expect_equal(sort(convert_2lattice_to_state(0,1/8,1,1,lat1,lat2)),lat_ref)
})

test_that("scaling shifting goes forward", {
lat1=c(1:4)/4-1/4
lat2=c(1:4)/4-1/4

lat_ref = c(1:8)/8-1/8
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(1/128,0,1,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,1/128,1,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,0,1.01,1,lat1,lat2)),)
expect_lt(max(convert_2lattice_to_state(0,0,1,1,lat1,lat2)),
          max(convert_2lattice_to_state(0,0,1,1.01,lat1,lat2)),)
})