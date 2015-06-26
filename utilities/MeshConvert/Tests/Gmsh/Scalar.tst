<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh box with scalar surface deformation and hex boundary layer splitting</description>
    <executable>MeshConvert</executable>
    <parameters>-m scalar:surf=1:scalar=exp\(-x*x-y*y\):nq=7 -m bl:nq=7:layers=7:r=2.5 -m jac:list Scalar.msh Scalar.xml:xml:test</parameters>
    <files>
        <file description="Input File">Scalar.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
