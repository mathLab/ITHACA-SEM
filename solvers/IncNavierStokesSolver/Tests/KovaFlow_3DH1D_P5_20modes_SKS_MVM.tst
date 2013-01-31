<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, P=5, 20 Fourier modes - Skew-Symmetric advection(MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_P5_20modes_SKS_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.4986e-06</value>
            <value variable="v" tolerance="1e-12">9.86701e-07</value>
            <value variable="w" tolerance="1e-12">9.32569e-07</value>
	    <value variable="p" tolerance="1e-12">2.08094e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">2.87017e-06</value>
            <value variable="v" tolerance="1e-12">2.25226e-06</value>
            <value variable="w" tolerance="1e-12">1.95419e-06</value>
	    <value variable="p" tolerance="1e-12">6.85396e-05</value>
        </metric>
    </metrics>
</test>