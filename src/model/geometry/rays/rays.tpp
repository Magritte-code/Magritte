///  Getter for the ray direction vector
///    @param[in]      o : number of point from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
////////////////////////////////////////////////////////////////////////
accel inline Vector3D Rays :: get_direction (
    const Size o,
    const Size r ) const
{
    Vector3D dir;

    if (parameters.adaptive_ray_tracing())
    {
        dir = m_direction(o, r);
    }
    else
    {
        dir = direction[r];
    }

    return dir;
}




///  Getter for the ray weight
///    @param[in]      o : number of point from which the ray originates
///    @param[in]      r : number of the ray along which we are looking
////////////////////////////////////////////////////////////////////////
accel inline Real Rays :: get_weight (
    const Size o,
    const Size r ) const
{
    Real wgt;

    if (parameters.adaptive_ray_tracing())
    {
        wgt = m_weight(o, r);
    }
    else
    {
        wgt = weight[r];
    }

    return wgt;
}
