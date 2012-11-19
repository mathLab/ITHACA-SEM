<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Laminar Channel Flow 3D homogeneous 1D, P=3, 20 Fourier modes (MVM)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_MVM.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">3.34833e-16</value>
            <value variable="v" tolerance="1e-6">1.21337e-16</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">2.64454e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.94289e-15</value>
            <value variable="v" tolerance="1e-6">4.75925e-16</value>
            <value variable="w" tolerance="1e-6">1.02571e-17</value>
            <value variable="p" tolerance="1e-6">1.4011e-13</value>
        </metric>
    </metrics>
</test>


