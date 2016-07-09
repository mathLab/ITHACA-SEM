<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=4-9, 20 Fourier modes with HDF5 input - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</file>
        <file description="Restart File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.10175e-05</value>
            <value variable="v" tolerance="1e-12">6.64294e-06</value>
            <value variable="w" tolerance="1e-12">5.74527e-06</value>
            <value variable="p" tolerance="1e-8">6.43048e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.31884e-05</value>
            <value variable="v" tolerance="1e-12">2.5058e-05</value>
            <value variable="w" tolerance="1e-12">1.70132e-05</value>
            <value variable="p" tolerance="1e-8">0.000138713</value>
        </metric>
    </metrics>
</test>
