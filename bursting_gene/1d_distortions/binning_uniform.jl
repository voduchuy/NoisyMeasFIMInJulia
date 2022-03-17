export Binning 

struct Binning <: DistortionModel1D
    width::Integer 
end

function conditional_probabilities(model::Binning, x::Integer, y::AbstractVector{<:Integer})
    out = zeros(Float64, length(y))
    @simd for i in 1:length(y)
        if (model.width*y[i] ≤ x < model.width*(y[i]+1))
             out[i] = 1.0
        end
    end
    out 
end

# def getConditionalProbabilities(self, x: int, y: np.ndarray) -> np.ndarray:
# return np.array((y*self.width <= x) & ((y+1)*self.width > x), dtype=np.double)