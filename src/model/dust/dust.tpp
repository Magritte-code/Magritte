///  Getter for the dust opacity
///    @param[in] p: point index
///    @param[in] freq: frequency
/////////////////////////////////
inline Real Dust :: get_opacity (const Size p, const Real freq) const
{
    Real opacity = density[p];

    const Size f = search (freqs, freq);

    if (f == 0)
    {
        opacity *= kappa[0];
    }
    else if (f == freqs.size()-1)
    {
        opacity *= kappa[freqs.size()-1];
    }
    else
    {
        const Real step = (freq - freqs[f-1]) / (freqs[f] - freqs[f-1]);

        opacity *= kappa[f-1] + (kappa[f] - kappa[f-1]) * step;
    }

    return opacity;
}
