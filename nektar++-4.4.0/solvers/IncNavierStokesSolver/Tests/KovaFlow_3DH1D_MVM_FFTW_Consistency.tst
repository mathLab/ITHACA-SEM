<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow 3D homogeneous 1D, Check MVM/FFTW consistency</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_3DH1D_MVM_FFTW_Consistency.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_3DH1D_MVM_FFTW_Consistency.xml</file>
        <file description="Restart file">KovaFlow_3DH1D_MVM_FFTW_Consistency.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.5498e-06</value>
            <value variable="v" tolerance="1e-08">9.5421e-07</value>
            <value variable="w" tolerance="1e-08">9.05058e-07</value>
	    <value variable="p" tolerance="1e-08">2.04453e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">3.14595e-06</value>
            <value variable="v" tolerance="1e-08">2.25562e-06</value>
            <value variable="w" tolerance="1e-08">1.81458e-06</value>
	    <value variable="p" tolerance="1e-08">6.87311e-05</value>
        </metric>
    </metrics>
</test>
