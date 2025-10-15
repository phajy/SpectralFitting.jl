using SpectralFitting, Test

grouping = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0]

items = [i for i in SpectralFitting.GroupingIterator(grouping)]
@test items == [(1, 1, 4), (2, 5, 7), (3, 8, 8), (4, 9, 10)]

items = collect(SpectralFitting.GroupingIterator(grouping))
@test items == [(1, 1, 4), (2, 5, 7), (3, 8, 8), (4, 9, 10)]

data = collect(range(0.0, 5.0, 10))
SpectralFitting.regroup!(data, grouping)

# should modify inplace
@test data ≈ [3.3333333333333335, 8.333333333333334, 3.888888888888889, 9.444444444444445]

# ensure having a 1 at the end doesn't break anything
grouping = [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1]

items = collect(SpectralFitting.GroupingIterator(grouping))
@test items == [(1, 1, 4), (2, 5, 7), (3, 8, 8), (4, 9, 11), (5, 12, 12)]

# Test signal-to-noise grouping
# Create a simple spectrum with known counts
# Using counts that will produce predictable S/N values
channels = collect(1:10)
quality = zeros(Int, 10)
grouping = ones(Int, 10)

# Create counts: [100, 100, 100, 25, 25, 25, 25, 9, 9, 9]
# S/N without background: sqrt(100)=10, sqrt(200)=14.14, sqrt(300)=17.32
# sqrt(25)=5, sqrt(50)=7.07, sqrt(75)=8.66, sqrt(100)=10
# sqrt(9)=3, sqrt(18)=4.24, sqrt(27)=5.20
data = Float64[100, 100, 100, 25, 25, 25, 25, 9, 9, 9]
errors = sqrt.(data)

spectrum = SpectralFitting.Spectrum{Float64}(
    channels,
    quality,
    grouping,
    data,
    u"counts",
    1.0,  # exposure_time
    1.0,  # background_scale
    1.0,  # area_scale
    SpectralFitting.ErrorStatistics.Poisson,
    errors,
    0.0,  # systematic_error
    "TEST",
    "TEST"
)

# Test grouping without background to S/N = 8
SpectralFitting.group_min_snr!(spectrum, 8.0)

# Expected behavior:
# Ch1: 100 counts, S/N=10 >= 8 → NEW_GRP (1)
# Ch2: 100 counts, S/N=10 >= 8 → NEW_GRP (1)
# Ch3: 100 counts, S/N=10 >= 8 → NEW_GRP (1)
# Ch4: 25 counts, S/N=5 < 8 → CONTINUE (0)
# Ch5: 25+25=50, S/N=7.07 < 8 → CONTINUE (0)
# Ch6: 50+25=75, S/N=8.66 >= 8 → NEW_GRP (1)
# Ch7-10: Continue accumulating
@test spectrum.grouping[1] == 1  # First channel starts a group
@test spectrum.grouping[2] == 1  # Second channel also meets threshold, starts new group
@test spectrum.grouping[3] == 1  # Third channel also meets threshold, starts new group
@test spectrum.grouping[4] == 0  # Fourth channel doesn't meet threshold, continues
@test spectrum.grouping[5] == 0  # Fifth channel accumulates
@test spectrum.grouping[6] == 1  # Sixth channel reaches threshold (25+25+25=75, S/N=8.66)
@test count(==(1), spectrum.grouping) >= 4  # At least 4 groups formed

# Test with minimum S/N = 15 (should require more counts per bin)
grouping2 = ones(Int, 10)
spectrum2 = SpectralFitting.Spectrum{Float64}(
    channels,
    quality,
    grouping2,
    copy(data),
    u"counts",
    1.0,
    1.0,
    1.0,
    SpectralFitting.ErrorStatistics.Poisson,
    copy(errors),
    0.0,
    "TEST",
    "TEST"
)

SpectralFitting.group_min_snr!(spectrum2, 15.0)

# S/N=15 requires at least 225 counts (15^2)
# Ch1-3: 100+100+100=300 gives S/N=17.32 >= 15 → grouping[3]=1, others=0
# Remaining channels don't accumulate enough
@test spectrum2.grouping[1] == 0  # Accumulating
@test spectrum2.grouping[2] == 0  # Still accumulating
@test spectrum2.grouping[3] == 1  # Threshold met (300 counts, S/N=17.32)
@test count(==(1), spectrum2.grouping) <= count(==(1), spectrum.grouping)  # Fewer groups with higher threshold

# Test with background
background_data = Float64[10, 10, 10, 5, 5, 5, 5, 2, 2, 2]
background = SpectralFitting.Spectrum{Float64}(
    channels,
    quality,
    ones(Int, 10),
    background_data,
    u"counts",
    1.0,
    1.0,
    1.0,
    SpectralFitting.ErrorStatistics.Poisson,
    sqrt.(background_data),
    0.0,
    "TEST",
    "TEST"
)

grouping3 = ones(Int, 10)
spectrum3 = SpectralFitting.Spectrum{Float64}(
    channels,
    quality,
    grouping3,
    copy(data),
    u"counts",
    1.0,
    1.0,
    1.0,
    SpectralFitting.ErrorStatistics.Poisson,
    copy(errors),
    0.0,
    "TEST",
    "TEST"
)

SpectralFitting.group_min_snr!(spectrum3, 5.0; background = background)

# With background, S/N calculation: (S - B*areanorm) / sqrt(S + B*areanorm^2)
# areanorm = 1.0 (same exposure and backscal)
# Ch1: (100 - 10*1) / sqrt(100 + 10*1^2) = 90 / sqrt(110) = 8.58 >= 5 → NEW_GRP (1)
# Ch2: (100 - 10*1) / sqrt(100 + 10*1^2) = 90 / sqrt(110) = 8.58 >= 5 → NEW_GRP (1)
# Ch3: Similar, S/N = 8.58 >= 5 → NEW_GRP (1)
@test spectrum3.grouping[1] == 1
@test spectrum3.grouping[2] == 1
@test spectrum3.grouping[3] == 1
@test count(==(1), spectrum3.grouping) >= 3
