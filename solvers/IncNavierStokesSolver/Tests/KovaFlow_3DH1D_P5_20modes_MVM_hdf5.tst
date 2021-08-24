<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=4-9, 20 Fourier modes in parallel with HDF5 input - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.xml</file>
        <file description="Restart File">KovaFlow_3DH1D_P5_20modes_MVM_hdf5.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.10128e-05</value>
            <value variable="v" tolerance="1e-12">6.64551e-06</value>
            <value variable="w" tolerance="1e-12">5.74741e-06</value>
            <value variable="p" tolerance="1e-8">6.45012e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.3134e-05</value>
            <value variable="v" tolerance="1e-12">2.51101e-05</value>
            <value variable="w" tolerance="1e-12">1.70084e-05</value>
            <value variable="p" tolerance="1e-8">0.000138221</value>
        </metric>
    </metrics>
</test>
