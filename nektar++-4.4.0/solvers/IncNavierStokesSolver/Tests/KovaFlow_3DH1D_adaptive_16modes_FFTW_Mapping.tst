<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, adaptive P, 16 Fourier modes, using explicit mapping</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-9">4.33081e-05</value>
            <value variable="v" tolerance="1e-9">1.16106e-05</value>
            <value variable="w" tolerance="1e-9">1.29757e-05</value>
	    <value variable="p" tolerance="1e-9">0.000260725</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">8.9668e-05</value>
            <value variable="v" tolerance="1e-9">2.47702e-05</value>
            <value variable="w" tolerance="1e-9">2.86236e-05</value>
	    <value variable="p" tolerance="1e-9">0.000805455</value>
        </metric>
    </metrics>
</test>
