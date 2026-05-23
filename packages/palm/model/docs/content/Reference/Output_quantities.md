# Output Quantities

---

## Vertical profile quantities
{{ include_palm_output_quantities('output_quantities', scopes=['all'], types=['vertical profile'], show_remarks=True) }}

## 3d-array quantities
{{ include_palm_output_quantities('output_quantities', scopes=['all'], types=['3d-array'], show_remarks=True) }}

## Masked array quantities
{{ include_palm_output_quantities('output_quantities', scopes=['all'], types=['masked array'], show_remarks=True) }}

## 2d-array quantities
{{ include_palm_output_quantities('output_quantities', scopes=['all'], types=['2d-array'], show_remarks=True) }}

## DET quantities
{{ include_palm_output_quantities('output_quantities', scopes=['det_model'], types=['3d-array', '2d-array', 'vertical profile'], show_remarks=True) }}

## SLUrb quantities
{{ include_palm_output_quantities('output_quantities', scopes=['slurb_model'], types=['3d-array', '2d-array'], show_remarks=True) }}

## Transport quantities
{{ include_palm_output_quantities('output_quantities', scopes=['traffic_model'], types=['3d-array'], show_remarks=True) }}

## UV quantities
{{ include_palm_output_quantities('output_quantities', scopes=['uv_radiation'], types=['2d-array'], show_remarks=True) }}